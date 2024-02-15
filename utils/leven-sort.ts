//@ts-nocheck
/**
 * leven-sort has an issue with being imported for some reason, since the module is somewhat small and stable
 *  I've just extracted the applicable code and put it here.
 */

/** leven
 * credit to https://github.com/sindresorhus/leven/blob/d3e23a590fedd726a3bcccc40739c34bfab18250/index.js
 * MIT License
 **/

const array = [];
const characterCodeCache = [];

function levenSrc(first, second) {
	if (first === second) {
		return 0;
	}

	const swap = first;

	// Swapping the strings if `a` is longer than `b` so we know which one is the
	// shortest & which one is the longest
	if (first.length > second.length) {
		first = second;
		second = swap;
	}

	let firstLength = first.length;
	let secondLength = second.length;

	// Performing suffix trimming:
	// We can linearly drop suffix common to both strings since they
	// don't increase distance at all
	// Note: `~-` is the bitwise way to perform a `- 1` operation
	while (firstLength > 0 && (first.charCodeAt(~-firstLength) === second.charCodeAt(~-secondLength))) {
		firstLength--;
		secondLength--;
	}

	// Performing prefix trimming
	// We can linearly drop prefix common to both strings since they
	// don't increase distance at all
	let start = 0;

	while (start < firstLength && (first.charCodeAt(start) === second.charCodeAt(start))) {
		start++;
	}

	firstLength -= start;
	secondLength -= start;

	if (firstLength === 0) {
		return secondLength;
	}

	let bCharacterCode;
	let result;
	let temporary;
	let temporary2;
	let index = 0;
	let index2 = 0;

	while (index < firstLength) {
		characterCodeCache[index] = first.charCodeAt(start + index);
		array[index] = ++index;
	}

	while (index2 < secondLength) {
		bCharacterCode = second.charCodeAt(start + index2);
		temporary = index2++;
		result = index2;

		for (index = 0; index < firstLength; index++) {
			temporary2 = bCharacterCode === characterCodeCache[index] ? temporary : temporary + 1;
			temporary = array[index];
			// eslint-disable-next-line no-multi-assign
			result = array[index] = temporary > result ? (temporary2 > result ? result + 1 : temporary2) : (temporary2 > temporary ? temporary + 1 : temporary2);
		}
	}

	return result;
}

/** leven-sort
 * credit to https://github.com/doesdev/leven-sort/blob/90497a37c9869d11a76a1dce010075db040ac327/index.js
 * MIT License
 **/

const bufChar = String.fromCharCode(2000)

export default (ary: string[], src1: string, key1?: string, src2?: string, key2?: string): string[] => {
  let max = 0

  ary.forEach(function (el) {
    if (!el) return
    if (!key1 && !key2 && el.length > max) max = el.length

    if (key1 instanceof Array) {
      return key1.forEach(function (k) {
        if (k && el[k] && el[k].length > max) max = el[k].length
      })
    }

    if (key1 && el[key1] && el[key1].length > max) max = el[key1].length
    if (key2 && el[key2] && el[key2].length > max) max = el[key2].length
  })

  const maximize = function (val) {
    if (!val || val.length === max) return val
    return val + Array(max - val.length).join(bufChar)
  }

  const maxDist = function (a, b, c) {
    const aLen = (a || '').length
    const bLen = (b || '').length

    return (aLen > bLen ? aLen : bLen) || 1
  }

  const leven = function (a, b) {
    if (a === b) return 0

    b = maximize(b)

    if (!a || !b) return maxDist(a, b)

    return levenSrc(a, b)
  }

  const levSort = function (src, a, b) {
    if (!a) return 1
    if (!b) return -1

    a = leven(src, a)
    b = leven(src, b)

    return a - b
  }

  const levMinInAry = function (array, src) {
    let min = 1000
    const len = array.length

    for (let counter = 0; counter < len; counter++) {
      const val = array[counter]

      if (val && val.length && val.length > 0) {
        const levScore = leven(src, array[counter])
        if (levScore < min) min = levScore
      }
    }

    return min
  }

  const sorted = ary.sort(function (a, b) {
    if (key1 instanceof Array) {
      const aLev = levMinInAry(key1.map(function (k) { return a[k] }), src1)
      const bLev = levMinInAry(key1.map(function (k) { return b[k] }), src1)

      return aLev - bLev
    }

    if (!key1 && !key2) return levSort(src1, a, b)

    if (!key2) return levSort(src1, a[key1], b[key1])

    const score = levSort(src1, a[key1], b[key1]) * 10

    return score + levSort(src2, a[key2], b[key2])
  })

  return sorted
}
