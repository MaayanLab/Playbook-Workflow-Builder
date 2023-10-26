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

export function mean_std(arr: number[]) {
  const mu = mean(arr)
  const std = Math.sqrt(sum(arr.map(v => Math.pow(v - mu, 2))) / arr.length)
  return { mu, std }
}
