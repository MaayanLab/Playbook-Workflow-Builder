import React from 'react'
import type { FileInputProps } from '@blueprintjs/core'
import * as array from '@/utils/array'
import dynamic from 'next/dynamic'

const Bp5FileInput = dynamic(() => import('@blueprintjs/core').then(({ FileInput }) => FileInput))

export default function FileInput(props: FileInputProps) {
  const [currentFilename, setCurrentFilename] = React.useState<string | undefined>(undefined)
  const props_ = {...props}
  if (!('draggable' in props)) props_.draggable = true
  props_.inputProps = 'inputProps' in props ? {...props.inputProps} : {}
  if (!('multiple' in props_.inputProps)) props_.inputProps.multiple = true
  const onChange = props_.inputProps.onChange
  props_.inputProps.onChange = (evt) => {
    if (!evt.currentTarget || !evt.currentTarget.files || evt.currentTarget.files.length === 0) {
      setCurrentFilename(undefined)
    } else {
      setCurrentFilename(
        array.arange(evt.currentTarget.files.length)
          .map(i => evt.currentTarget.files && evt.currentTarget.files.item(i))
          .filter((file): file is File => file !== null)
          .map(file => file.name)
          .join(', ')
      )
    }
    if (onChange) onChange(evt)
  }
  if (!('text' in props)) props_.text = currentFilename
  if (!('hasSelection' in props)) props_.hasSelection = !!currentFilename
  return <Bp5FileInput {...props_} />
}
