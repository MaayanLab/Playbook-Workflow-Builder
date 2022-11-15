import React from 'react'
import * as array from '@/utils/array'
import { Table2 as Table, Column, RowHeaderCell2, Cell } from '@blueprintjs/table'

interface MatrixProps {
  /**
   * Labels on the columns of the matrix
   */
  columns: Array<string>,
  /**
   * Labels on the row of the matrix
   */
  index: Array<string>,
  /**
   * Actual values of the matrix (array of arrays)
   */
  values: Array<Array<unknown>>,
  /**
   * [optional]: Where to place ellipses, when showing a truncated matrix
   * of the form: [row-elipse-index | null, column-elipse-index | null]
   */
  ellipses?: [number | undefined, number | undefined]
}

export default function Matrix(props: MatrixProps) {
  return (
    <Table
      numRows={props.values.length + (props.ellipses && props.ellipses[0] ? 1 : 0)}
      enableRowResizing
      enableColumnResizing
      rowHeaderCellRenderer={(row) => {
        let row_ellipse = false
        if (props.ellipses && props.ellipses[0]) {
          if (row === props.ellipses[0]) row_ellipse = true
          else if (row > props.ellipses[0]) row = row - 1
        }
        return (
          <RowHeaderCell2
            key={row_ellipse ? '-ellipse' : (row+'')}
            name={row_ellipse ? '...' : (props.index[row]+'')}
          />
        )
      }}
    >
      {array.arange(props.columns.length + (props.ellipses && props.ellipses[1] ? 1 : 0)).map((column) => {
        let column_ellipse = false
        if (props.ellipses && props.ellipses[1]) {
          if (column === props.ellipses[1]) column_ellipse = true
          else if (column > props.ellipses[1]) column = column - 1
        }
        return (
          <Column
            key={`${column}${column_ellipse ? '-ellipse' : ''}`}
            name={column_ellipse ? '...' : props.columns[column]}
            cellRenderer={(row, _) => {
              let row_ellipse = column_ellipse
              if (props.ellipses && props.ellipses[0]) {
                if (row === props.ellipses[0]) row_ellipse = true
                else if (row > props.ellipses[0]) row = row - 1
              }
              return (
                <Cell key={row_ellipse ? '-ellipse' : (row+'')}>
                  {row_ellipse ? '...' : (props.values[row][column]+'')}
                </Cell>
              )
            }}
          />
        )
      })}
    </Table>
  )
}
