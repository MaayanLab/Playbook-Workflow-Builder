/**
 * A prefix tree for identification of the correct
 *  split given a list of path prefixes.
 * 
 * e.g.
 * const t = create_prefix_tree_from_paths([
 *  'a/b/c',
 *  'a/d',
 * ])
 * search_prefix_tree(t, 'a/d/e') == {
 *   prefix: 'a/d',
 *   path: 'e'
 * }
 */

type TreeT<T> = Record<string, T | string>
interface Tree extends TreeT<Tree> {}

export function create_prefix_tree_from_paths(paths: string[]) {
  const tree: Tree = {}
  for (const path of paths) {
    const path_split = path.split('/')
    let tmp = tree
    for (let i = 0; i < path_split.length-1; i++) {
      if (!(path_split[i] in tmp)) {
        tmp[path_split[i]] = {}
      }
      tmp = tmp[path_split[i]] as Tree
    }
    tmp[path_split[path_split.length-1]] = path
  }
  return tree
}

export function search_prefix_tree(tree: Tree, path: string) {
  const path_split = path.split('/')
  let prefix: string | undefined, prefix_path: string | undefined
  let tmp: Tree = tree
  for (let i = 0; i <= path_split.length; i++) {
    const part = path_split[i]
    if (part in tmp) {
      tmp = tmp[part] as Tree
      if (typeof tmp === 'string') {
        prefix = tmp
        prefix_path = path_split.slice(i+1).join('/')
        break
      }
    } else break
  }
  return { prefix, path: prefix_path }
}
