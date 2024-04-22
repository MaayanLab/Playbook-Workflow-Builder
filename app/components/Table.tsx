import { Table2 as Bp5Table, Column as Bp5Column, Cell as Bp5Cell, RowHeaderCell2 as Bp5RowHeaderCell, EditableCell2 as Bp5EditableCell, Table2Props as Bp5Table2Props } from '@blueprintjs/table'
import * as dict from '@/utils/dict'
import dynamic from 'next/dynamic'

const Bp5Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))
const Bp5Menu = dynamic(() => import('@blueprintjs/core').then(({ Menu }) => Menu))
const Bp5MenuItem = dynamic(() => import('@blueprintjs/core').then(({ MenuItem }) => MenuItem))
const Bp5Popover = dynamic(() => import('@blueprintjs/core').then(({ Popover }) => Popover))

export const Cell = Bp5Cell
export const Column = Bp5Column
export const RowHeaderCell = Bp5RowHeaderCell
export const EditableCell = Bp5EditableCell
export type Table2Props = Bp5Table2Props & {
  height?: number
  downloads?: Record<string, () => void>,
  shape?: Array<number>
}

export function Table({ children, height, shape: shape_, downloads, ...props }: Table2Props) {
  const shape = shape_ ? shape_ : Array.isArray(children) ? [props.numRows, children.length] : [props.numRows]
  return (
    <>
      <div style={{ maxHeight: height, overflow: 'hidden' }}>
        <Bp5Table {...props}>{children}</Bp5Table>
      </div>
      <div className="bp5-navbar-group">
        <span className="prose max-w-none">Shape: ({shape.map((dim, i) => <span key={i}>{i>0?', ':''}{dim}</span>)})</span>
        <span className="bp5-navbar-divider"></span>
        {downloads ?
          dict.keys(downloads).length === 1 ?
            <Bp5Button text="Download" icon="download" minimal onClick={() => {
              if (!downloads) return
              const [download] = dict.values(downloads)
              download()
            }} />
          : <Bp5Popover
              content={
                <Bp5Menu>
                  {dict.items(downloads).map(({ key, value }) =>
                    <Bp5MenuItem key={key} icon="document" text={key} onClick={() => {value()}} />
                  )}
                </Bp5Menu>
              }
              placement="bottom"
            >
              <Bp5Button text="Download" icon="download" minimal />
            </Bp5Popover>
          : null}
      </div>
    </>
  )
}
