import React from 'react'
import Masonry from 'react-masonry-css'
import tsvector from "@/utils/tsvector"
import dynamic from "next/dynamic"
import * as array from '@/utils/array'
// import { useQsState } from "@/components/qsstate"

const FormGroup = dynamic(() => import('@blueprintjs/core').then(({ FormGroup }) => FormGroup))
const InputGroup = dynamic(() => import('@blueprintjs/core').then(({ InputGroup }) => InputGroup))

type KVCounts = { [key: string]: { [val: string]: number } }

export default function Catalog<T extends { meta?: { pagerank?: number, tags?: Record<string, Record<string, number>> } }>({ items, serialize, children }: {
  items: Array<T>,
  serialize: (item: T) => string,
  children: (child: T) => React.ReactElement,
}) {
  const [search, setSearch] = React.useState('')
  const [filters, setFilters] = React.useState({} as Record<string, Record<string, number>>)
  const _search = tsvector(search)
  const _group_values: KVCounts = {}
  for (const k in items) {
    const item = items[k]
    const item_meta = item.meta||{}
    if (!('tags' in item_meta)) continue
    for (const group in item_meta.tags) {
      if (!(group in _group_values)) {
        _group_values[group] = {}
      }
      for (const value in item_meta.tags[group]) {
        if (!(value in _group_values[group])) {
          _group_values[group][value] = 0
        }
        _group_values[group][value] += 1
      }
    }
  }
  // filtered
  const _items_filtered = Object.values(items)
    .filter((item) => {
      const item_meta = item.meta || {}
      const item_meta_tags = item_meta.tags || {}
      return array.all(
        Object.keys(filters)
          .map(group => array.any(
            Object.keys(item_meta_tags[group]||{})
              .map(value => (filters[group]||{})[value] === 1)
          ))
      )
    })
    .filter((item) => {
      if (_search.size === 0) return true
      const _fts = tsvector(serialize(item))
      return _fts.intersect(_search).size !== 0
    })
  // sort
  _items_filtered.sort((a, b) => {
    const a_meta = a.meta || {}
    const b_meta = b.meta || {}
    return (b_meta.pagerank||0) - (a_meta.pagerank||0)
  })

  const _group_values_filtered: KVCounts = {}
  for (const k in _items_filtered) {
    const item = _items_filtered[k]
    const item_meta = item.meta||{}
    const item_meta_tags = item_meta.tags || {}
    for (const group in item_meta_tags) {
      if (!(group in _group_values_filtered)) {
        _group_values_filtered[group] = {}
      }
      for (const value in item_meta_tags[group]) {
        if (!(value in _group_values_filtered[group])) {
          _group_values_filtered[group][value] = 0
        }
        _group_values_filtered[group][value] += 1
      }
    }
  }
  return (
    <div className="flex-grow flex flex-col">
      <FormGroup label="Filter">
        <InputGroup
          leftIcon="filter"
          placeholder="Filter string..."
          value={search}
          onChange={evt => setSearch(evt.currentTarget.value)}
        />
      </FormGroup>
      <div className="flex-grow flex flex-row">
        <div className="sm:col-span-12 md:col-span-4 lg:col-span-3">
          {Object.keys(_group_values)
            .filter(group => Object.keys(_group_values[group]).length > 1)
            .map(group => (
              <fieldset key={group} className="pr-2">
                <legend>{group}</legend>
                <div className="ml-2 mb-4">
                  {Object.keys(_group_values[group])
                    .map(value => (
                      <label key={value} className="bp4-control bp4-switch">
                        <input
                          type="checkbox"
                          checked={(filters[group] || {})[value] === 1}
                          onChange={evt => {
                            const _filters = {...filters}
                            if (!evt.currentTarget.checked) {
                              if (!(group in _filters)) _filters[group] = {}
                              else _filters[group] = {...filters[group]}
                              if (value in _filters[group]) {
                                delete _filters[group][value]
                              }
                              if (Object.keys(_filters[group]).length === 0) {
                                delete _filters[group]
                              }
                            } else {
                              if (!(group in _filters)) _filters[group] = {}
                              else _filters[group] = {...filters[group]}
                              if (!(value in _filters[group])) {
                                _filters[group][value] = 1
                              }
                            }
                            setFilters(_filters)
                          }}
                        />
                        <span className="bp4-control-indicator"></span>
                        {value} <span className="whitespace-nowrap">[{
                        (_group_values_filtered[group]||{})[value] === _group_values[group][value] ? (
                          _group_values[group][value]
                        ) : (
                          `${(_group_values_filtered[group]||{})[value]||0} / ${_group_values[group][value]}`
                        )} card{_group_values[group][value] > 1 ? 's' : ''}]</span>
                      </label>
                  ))}
                </div>
              </fieldset>
            )
          )}
        </div>
        <div className="flex-grow">
          <Masonry
            breakpointCols={{
              // note these breakpoints match the bootstrap breakpoints
              //  they should only be changed along with the columnClassName spec
              default: 4,
              1400: 4,
              1200: 3,
              992: 2,
              768: 1,
            }}
            className="flex flex-row gap-2"
            columnClassName="flex-grow flex flex-col gap-2"
          >
            {_items_filtered.map(item => children(item))}
          </Masonry>
        </div>
      </div>
    </div>
  )
}
