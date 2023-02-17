import { Button as Bp4Button, Menu as Bp4Menu, MenuDivider as Bp4MenuDivider, MenuItem as Bp4MenuItem } from '@blueprintjs/core'
import { Table2 as Bp4Table, Column as Bp4Column, Cell as Bp4Cell, Table2Props as Bp4Table2Props, RowHeaderCell2 as Bp4RowHeaderCell, EditableCell2 as Bp4EditableCell } from '@blueprintjs/table'
import { Popover2 as Bp4Popover } from '@blueprintjs/popover2'
import * as dict from '@/utils/dict'

export const Cell = Bp4Cell
export const Column = Bp4Column
export const RowHeaderCell = Bp4RowHeaderCell
export const EditableCell = Bp4EditableCell
export type Table2Props = Bp4Table2Props & {
  height?: number
  downloads?: Record<string, () => void>,
  shape?: Array<number>
}

export function Table({ children, height, shape: shape_, downloads, ...props }: Table2Props) {
  const shape = shape_ ? shape_ : Array.isArray(children) ? [props.numRows, children.length] : [props.numRows]
  return (
    <>
      <div style={{ height }}>
        <Bp4Table {...props}>{children}</Bp4Table>
      </div>
      <div className="bp4-navbar-group">
        <span>Shape: ({shape.map((dim, i) => <span key={i}>{i>0?', ':''}{dim}</span>)})</span>
        <Bp4MenuDivider />
        {downloads ?
          dict.keys(downloads).length === 1 ?
            <Bp4Button text="Download" icon="download" minimal onClick={() => {
              if (!downloads) return
              const [download] = dict.values(downloads)
              download()
            }} />
          : <Bp4Popover
              content={
                <Bp4Menu>
                  {dict.items(downloads).map(({ key, value }) =>
                    <Bp4MenuItem key={key} icon="document" text={key} onClick={() => {value()}} />
                  )}
                </Bp4Menu>
              }
              placement="bottom"
            >
              <Bp4Button text="Download" icon="download" minimal />
            </Bp4Popover>
          : null}
      </div>
    </>
  )
}
