/*
 * Project component modules are deliberately authored as cohesive MetaNode DSL
 * definitions: their React view/prompt and server resolver live together.  A
 * browser build must not retain the resolver closure (or imports used only by
 * that closure), however.
 *
 * This webpack loader projects those modules for the browser without changing
 * their source files.  Server builds never run this loader.
 */
const ts = require('typescript')

const CLIENT_RESOLVE_ERROR =
  'MetaNode resolvers execute on the server; use the resolution API from the browser'

function isResolveCall(node) {
  return ts.isCallExpression(node)
    && ts.isPropertyAccessExpression(node.expression)
    && node.expression.name.text === 'resolve'
    && node.arguments.length === 1
    && (
      ts.isArrowFunction(node.arguments[0])
      || ts.isFunctionExpression(node.arguments[0])
      || ts.isPropertyAccessExpression(node.arguments[0])
      || ts.isIdentifier(node.arguments[0])
    )
}

function clientResolveStub() {
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
          factory.createStringLiteral(CLIENT_RESOLVE_ERROR),
        ]),
      ),
    ], true),
  )
}

function projectResolvers(sourceFile) {
  const transformer = context => {
    const visit = node => {
      if (isResolveCall(node)) {
        return ts.factory.updateCallExpression(
          node,
          ts.visitNode(node.expression, visit),
          node.typeArguments,
          [clientResolveStub()],
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

function usedIdentifiers(sourceFile) {
  const used = new Set()
  const visit = node => {
    // Imports declare bindings; counting them would prevent their removal.
    if (ts.isImportDeclaration(node) || ts.isImportEqualsDeclaration(node)) return
    if (ts.isIdentifier(node)) used.add(node.text)
    ts.forEachChild(node, visit)
  }
  ts.forEachChild(sourceFile, visit)
  return used
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
    const kind = this.resourcePath.endsWith('.tsx') ? ts.ScriptKind.TSX : ts.ScriptKind.TS
    const sourceFile = ts.createSourceFile(this.resourcePath, source, ts.ScriptTarget.Latest, true, kind)
    const projected = removeUnusedImports(projectResolvers(sourceFile))
    const output = ts.createPrinter({ newLine: ts.NewLineKind.LineFeed }).printFile(projected)
    callback(null, output, inputSourceMap)
  } catch (error) {
    callback(error)
  }
}

module.exports.CLIENT_RESOLVE_ERROR = CLIENT_RESOLVE_ERROR
