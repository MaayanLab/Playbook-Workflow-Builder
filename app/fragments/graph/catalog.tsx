import React from 'react'
import Masonry from 'react-masonry-css'
import { tsvector, tsvector_intersect } from "@/utils/tsvector"
import dynamic from "next/dynamic"
import * as array from '@/utils/array'
import * as dict from '@/utils/dict'

const FormGroup = dynamic(() => import('@blueprintjs/core').then(({ FormGroup }) => FormGroup))
const InputGroup = dynamic(() => import('@blueprintjs/core').then(({ InputGroup }) => InputGroup))

type KVCounts = { [key: string]: { [val: string]: number } }

export default function Catalog<T extends { spec: string, meta?: { pagerank?: number, tags?: Record<string, Record<string, number>> } }>({ items, weights, serialize, children }: {
  items: Array<T>,
  weights: Record<string, number>,
  serialize: (item: T) => string,
  children: (child: T) => React.ReactElement,
}) {
  const [search, setSearch] = React.useState('')
  const [filters, setFilters] = React.useState({} as Record<string, Record<string, number>>)
  const search_ts = React.useMemo(() => tsvector(search), [search])
  const { group_values, item_search_ts, pagerank_max, weight_max } = React.useMemo(() => {
    const group_values: KVCounts = {}
    const item_search_ts: Record<string, Set<string>> = {}
    let pagerank_max = 1
    let weight_max = 1
    for (const k in items) {
      const item = items[k]
      const item_meta = item.meta||{}
      item_search_ts[k] = tsvector(serialize(item))
      pagerank_max = Math.max(item_meta.pagerank||0, pagerank_max)
      weight_max = Math.max(weights[item.spec]||0, weight_max)
      if (!('tags' in item_meta)) continue
      for (const group in item_meta.tags) {
        if (!(group in group_values)) {
          group_values[group] = {}
        }
        for (const value in item_meta.tags[group]) {
          if (!(value in group_values[group])) {
            group_values[group][value] = 0
          }
          group_values[group][value] += 1
        }
      }
    }
    return { group_values, item_search_ts, pagerank_max, weight_max }
  }, [items])
  const items_filtered = React.useMemo(() =>
    array.arange(items.length)
      .filter((k) => {
        const item = items[k]
        const item_meta = item.meta || {}
        const item_meta_tags = item_meta.tags || {}
        return array.all(
          dict.keys(filters)
            .map(group => array.any(
              dict.keys(item_meta_tags[group]||{})
                .map(value => (filters[group]||{})[value] === 1)
            ))
        )
      }),
    [items, filters]
  )
  const items_filtered_searched = React.useMemo(() => {
    let items_filtered_search_score_max = 1
    const items_filtered_search_score: Record<number, number> = {}
    const items_filtered_searched =
      !search ? items_filtered
      : items_filtered
        .filter((k) => {
          const _fts = item_search_ts[k]
          const score = tsvector_intersect(_fts, search_ts).size
          items_filtered_search_score[k] = score
          items_filtered_search_score_max = Math.max(score, items_filtered_search_score_max)
          return score !== 0
        })
    items_filtered_searched.sort((a, b) => {
      const a_meta = items[a].meta || {}
      const a_pagerank = (a_meta.pagerank||0.)/pagerank_max
      const a_weight = (weights[items[a].spec]||0.)/weight_max
      const b_meta = items[b].meta || {}
      const b_pagerank = (b_meta.pagerank||0.)/pagerank_max
      const b_weight = (weights[items[b].spec]||0.)/weight_max
      return (
        (b_pagerank
          + 2 * b_weight
          + 2 * ((items_filtered_search_score[b]||0.)/items_filtered_search_score_max))
        -
        (a_pagerank
          + 2 * a_weight
          + 2 * ((items_filtered_search_score[a]||0.)/items_filtered_search_score_max))
      )
    })
    return items_filtered_searched.map(k => items[k])
  },
    [items_filtered, search_ts, item_search_ts, pagerank_max]
  )
  const group_values_filtered = React.useMemo(() => {
    const group_values_filtered: KVCounts = {}
    for (const k in items_filtered_searched) {
      const item = items_filtered_searched[k]
      const item_meta = item.meta||{}
      const item_meta_tags = item_meta.tags || {}
      for (const group in item_meta_tags) {
        if (!(group in group_values_filtered)) {
          group_values_filtered[group] = {}
        }
        for (const value in item_meta_tags[group]) {
          if (!(value in group_values_filtered[group])) {
            group_values_filtered[group][value] = 0
          }
          group_values_filtered[group][value] += 1
        }
      }
    }
    return group_values_filtered
  }, [items_filtered_searched])
  return (
    <div className="flex-grow flex flex-col mx-4">
      <FormGroup label={<span className="prose">Filter</span>}>
        <InputGroup
          leftIcon="filter"
          placeholder="Filter string..."
          value={search}
          onChange={evt => setSearch(evt.currentTarget.value)}
        />
      </FormGroup>
      <div className="flex-grow flex flex-row flex-wrap sm:flex-nowrap">
        <div className="flex-none">
          {dict.keys(group_values)
            .filter(group => dict.keys(group_values[group]).length > 1)
            .map(group => (
              <fieldset key={group} className="pr-2">
                <legend>{group}</legend>
                <div className="ml-2 mb-4">
                  {dict.keys(group_values[group])
                    .map(value => (
                      <label key={value} className="bp5-control bp5-switch prose prose-sm whitespace-nowrap">
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
                              if (dict.isEmpty(_filters[group])) {
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
                        <span className="bp5-control-indicator"></span>
                        {value} ({
                        (group_values_filtered[group]||{})[value] === group_values[group][value] ? (
                          group_values[group][value]
                        ) : (
                          `${(group_values_filtered[group]||{})[value]||0} / ${group_values[group][value]}`
                        )})
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
              1279: 3,
              1023: 2,
              767: 1,
            }}
            className="flex flex-row gap-2"
            columnClassName="flex-grow flex flex-col gap-2"
          >
            {items_filtered_searched.map(item => children(item))}
          </Masonry>
        </div>
      </div>
    </div>
  )
}
