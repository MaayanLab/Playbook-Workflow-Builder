import React from 'react'
import ReactSupervenn from 'react-supervenn'
import styles from 'react-supervenn/dist/lib/index.module.css'

export default function Supervenn(props: Omit<React.ComponentProps<typeof ReactSupervenn>, 'style'>) {
  return <ReactSupervenn style={styles} {...props} />
}
