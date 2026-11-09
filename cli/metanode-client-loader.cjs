/*
 * Project component modules are deliberately authored as cohesive MetaNode DSL
 * definitions: their React view/prompt and server resolver live together.  A
 * browser build must not retain the resolver closure (or imports used only by
 * that closure), however.
 *
 * This (webpack/Turbopack) loader projects those modules per target without
 * changing their source files: the browser build drops resolvers, the server build
 * drops views & prompts (along with lazy `next/dynamic` components only they used,
 * `ssr: false` is not allowed in the server graph).
 */
const ts = require('typescript')

const PROJECTION_ERRORS = {
  resolve: 'MetaNode resolvers execute on the server; use the resolution API from the browser',
  view: 'MetaNode views render in the browser and are unavailable on the server',
  prompt: 'MetaNode prompts render in the browser and are unavailable on the server',
}

function isProjectedCall(node, methods) {
  return ts.isCallExpression(node)
    && ts.isPropertyAccessExpression(node.expression)
    && methods.has(node.expression.name.text)
    && node.arguments.length === 1
    && (
      ts.isArrowFunction(node.arguments[0])
      || ts.isFunctionExpression(node.arguments[0])
      || ts.isPropertyAccessExpression(node.arguments[0])
      || ts.isIdentifier(node.arguments[0])
    )
}

function projectedStub(method) {
  const factory = ts.factory
  return factory.createArrowFunction(
    [factory.createModifier(ts.SyntaxKind.AsyncKeyword)],
    undefined,
    [],
    undefined,
    factory.createToken(ts.SyntaxKind.EqualsGreaterThanToken),
    factory.createBlock([
      factory.createThrowStatement(
        factory.createNewExpression(factory.createIdentifier('Error'), undefined, [
          factory.createStringLiteral(PROJECTION_ERRORS[method]),
        ]),
      ),
    ], true),
  )
}

function projectMethods(sourceFile, methods, target) {
  const transformer = context => {
    const visit = node => {
      if (
        target === 'server'
        && ts.isCallExpression(node)
        && ts.isPropertyAccessExpression(node.expression)
        && node.expression.name.text === 'meta'
        && node.arguments.length === 1
        && ts.isObjectLiteralExpression(node.arguments[0])
      ) {
        const metadata = node.arguments[0]
        const properties = metadata.properties.filter(property => {
          if (!ts.isPropertyAssignment(property) && !ts.isShorthandPropertyAssignment(property)) return true
          const name = property.name
          return !(ts.isIdentifier(name) || ts.isStringLiteral(name)) || (name.text !== 'icon' && name.text !== 'logo')
        })
        return ts.factory.updateCallExpression(
          node,
          ts.visitNode(node.expression, visit),
          node.typeArguments,
          [ts.factory.updateObjectLiteralExpression(metadata, properties.map(property => ts.visitNode(property, visit)))],
        )
      }
      if (isProjectedCall(node, methods)) {
        const method = node.expression.name.text
        return ts.factory.updateCallExpression(
          node,
          ts.visitNode(node.expression, visit),
          node.typeArguments,
          [projectedStub(method)],
        )
      }
      return ts.visitEachChild(node, visit, context)
    }
    return root => ts.visitNode(root, visit)
  }
  const result = ts.transform(sourceFile, [transformer])
  const transformed = result.transformed[0]
  result.dispose()
  return transformed
}

function usedIdentifiers(sourceFile, { skipLazyComponents = false } = {}) {
  const used = new Set()
  const visit = node => {
    // a lazy component's initializer can mention its own name (e.g. `.then(({ X }) => X)`)
    if (skipLazyComponents && isLazyComponentDeclaration(node)) return
    // Imports declare bindings; counting them would prevent their removal.
    if (ts.isImportDeclaration(node) || ts.isImportEqualsDeclaration(node)) return
    // likewise a variable's own name is a binding, not a use of it
    if (ts.isVariableDeclaration(node) && ts.isIdentifier(node.name)) {
      if (node.type) visit(node.type)
      if (node.initializer) visit(node.initializer)
      return
    }
    if (ts.isIdentifier(node)) used.add(node.text)
    ts.forEachChild(node, visit)
  }
  ts.forEachChild(sourceFile, visit)
  return used
}

/** `const X = dynamic(...)` (not exported) */
function isLazyComponentDeclaration(statement) {
  return ts.isVariableStatement(statement)
    && !statement.modifiers?.some(modifier => modifier.kind === ts.SyntaxKind.ExportKeyword)
    && statement.declarationList.declarations.every(declaration =>
      ts.isIdentifier(declaration.name)
      && declaration.initializer
      && ts.isCallExpression(declaration.initializer)
      && ts.isIdentifier(declaration.initializer.expression)
      && declaration.initializer.expression.text === 'dynamic'
    )
}

/**
 * Remove lazy components which are no longer referenced (they were only used by
 *  a projected view/prompt), creating a dynamic import is free of side effects.
 */
function removeUnusedLazyComponents(sourceFile) {
  const used = usedIdentifiers(sourceFile, { skipLazyComponents: true })
  const statements = sourceFile.statements.filter(statement =>
    !isLazyComponentDeclaration(statement)
    || statement.declarationList.declarations.some(declaration => used.has(declaration.name.text))
  )
  return ts.factory.updateSourceFile(sourceFile, statements)
}

function removeUnusedImports(sourceFile) {
  const used = usedIdentifiers(sourceFile)
  const statements = sourceFile.statements.flatMap(statement => {
    if (!ts.isImportDeclaration(statement) || !statement.importClause) return [statement]

    const clause = statement.importClause
    const name = clause.name && used.has(clause.name.text) ? clause.name : undefined
    let namedBindings
    if (clause.namedBindings && ts.isNamespaceImport(clause.namedBindings)) {
      if (used.has(clause.namedBindings.name.text)) namedBindings = clause.namedBindings
    } else if (clause.namedBindings && ts.isNamedImports(clause.namedBindings)) {
      const elements = clause.namedBindings.elements.filter(element => used.has(element.name.text))
      if (elements.length) namedBindings = ts.factory.updateNamedImports(clause.namedBindings, elements)
    }

    // Preserve explicit side-effect imports, but an import whose bindings all
    // disappeared was only supporting a stripped resolver and can be removed.
    if (!name && !namedBindings) return []
    return [ts.factory.updateImportDeclaration(
      statement,
      statement.modifiers,
      ts.factory.updateImportClause(clause, clause.isTypeOnly, name, namedBindings),
      statement.moduleSpecifier,
      statement.attributes,
    )]
  })
  return ts.factory.updateSourceFile(sourceFile, statements)
}

module.exports = function metanodeClientLoader(source, inputSourceMap) {
  const callback = this.async()
  try {
    if (!source.includes('MetaNode(')) return callback(null, source, inputSourceMap)
    const target = this.getOptions?.().target ?? 'client'
    const methods = new Set(target === 'server' ? ['view', 'prompt'] : ['resolve'])
    const kind = this.resourcePath.endsWith('.tsx') ? ts.ScriptKind.TSX : ts.ScriptKind.TS
    const sourceFile = ts.createSourceFile(this.resourcePath, source, ts.ScriptTarget.Latest, true, kind)
    const projected = removeUnusedImports(removeUnusedLazyComponents(projectMethods(sourceFile, methods, target)))
    const output = ts.createPrinter({ newLine: ts.NewLineKind.LineFeed }).printFile(projected)
    callback(null, output, inputSourceMap)
  } catch (error) {
    callback(error)
  }
}

module.exports.PROJECTION_ERRORS = PROJECTION_ERRORS
