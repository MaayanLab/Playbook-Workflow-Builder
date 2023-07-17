export function sum(arr: number[]) {
  let acc = 0
  for (let el of arr) {
    acc += el
  }
  return acc
}

export function mean(arr: number[]) {
  return sum(arr) / arr.length
}

export function absmax([first, ...arr]: number[]) {
  let max = first
  for (let el of arr) {
    if (Math.abs(el) > max) {
      max = el
    }
  }
  return max
}
