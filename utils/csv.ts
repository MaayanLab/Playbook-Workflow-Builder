export default function read_csv<C extends string>(content: string, delimiter = ',') {
  const content_split = content.split(/\n/g)
  const header = content_split.shift()
  const columns = (header ? header.split(delimiter) : []) as C[]
  return {
    columns,
    values: content_split
      .map(row => row.split(delimiter))
      .map(row => {
        const record: { [k: string]: string } = {}
        columns.forEach((col, ind) => {
          record[col] = row[ind]
        })
        return record as Record<C, string>
      }),
  }
}