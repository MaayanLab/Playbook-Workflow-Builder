import { Table2 as Bp4Table, Column as Bp4Column, Cell as Bp4Cell, RowHeaderCell2 as Bp4RowHeaderCell, EditableCell2 as Bp4EditableCell, Table2Props as Bp4Table2Props } from '@blueprintjs/table'
import * as dict from '@/utils/dict'
import dynamic from 'next/dynamic'

const Bp4Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))
const Bp4Menu = dynamic(() => import('@blueprintjs/core').then(({ Menu }) => Menu))
const Bp4MenuDivider = dynamic(() => import('@blueprintjs/core').then(({ MenuDivider }) => MenuDivider))
const Bp4MenuItem = dynamic(() => import('@blueprintjs/core').then(({ MenuItem }) => MenuItem))
const Bp4Popover = dynamic(() => import('@blueprintjs/popover2').then(({ Popover2 }) => Popover2))

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
