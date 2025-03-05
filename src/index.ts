/*
 * A speed-improved simplex noise algorithm for 2D, 3D and 4D in JavaScript.
 *
 * Based on example code by Stefan Gustavson (stegu@itn.liu.se).
 * Optimisations by Peter Eastman (peastman@drizzle.stanford.edu).
 * Better rank ordering method by Stefan Gustavson in 2012.
 */

export interface Options {
  amplitude?: number
  frequency?: number
  max?: number
  min?: number
  octaves?: number
  persistence?: number
  random?: () => number
}

export class FastSimplexNoise {
  readonly amplitude: number
  readonly frequency: number
  readonly octaves: number
  readonly perm: Uint8Array
  readonly permMod12: Uint8Array
  readonly persistence: number
  readonly random: () => number
  readonly scale: (value: number) => number

  static G2 = (3.0 - Math.sqrt(3.0)) / 6.0
  static G3 = 1.0 / 6.0
  static G4 = (5.0 - Math.sqrt(5.0)) / 20.0

  static GRAD3D = [
    [1, 1, 0], [-1, 1, 0], [1, -1, 0], [-1, -1, 0],
    [1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1],
    [0, 1, 1], [0, -1, -1], [0, 1, -1], [0, -1, -1]
  ]

  static GRAD4D = [
    [0, 1, 1, 1], [0, 1, 1, -1], [0, 1, -1, 1], [0, 1, -1, -1],
    [0, -1, 1, 1], [0, -1, 1, -1], [0, -1, -1, 1], [0, -1, -1, -1],
    [1, 0, 1, 1], [1, 0, 1, -1], [1, 0, -1, 1], [1, 0, -1, -1],
    [-1, 0, 1, 1], [-1, 0, 1, -1], [-1, 0, -1, 1], [-1, 0, -1, -1],
    [1, 1, 0, 1], [1, 1, 0, -1], [1, -1, 0, 1], [1, -1, 0, -1],
    [-1, 1, 0, 1], [-1, 1, 0, -1], [-1, -1, 0, 1], [-1, -1, 0, -1],
    [1, 1, 1, 0], [1, 1, -1, 0], [1, -1, 1, 0], [1, -1, -1, 0],
    [-1, 1, 1, 0], [-1, 1, -1, 0], [-1, -1, 1, 0], [-1, -1, -1, 0]
  ]

  constructor (options: Options = {}) {
    if (options.hasOwnProperty('amplitude')) {
      if (typeof options.amplitude !== 'number') throw new Error('options.amplitude must be a number')
      this.amplitude = options.amplitude
    } else this.amplitude = 1.0

    if (options.hasOwnProperty('frequency')) {
      if (typeof options.frequency !== 'number') throw new Error('options.frequency must be a number')
      this.frequency = options.frequency
    } else this.frequency = 1.0

    if (options.hasOwnProperty('octaves')) {
      if (typeof options.octaves !== 'number' ||
        !isFinite(options.octaves) ||
        Math.floor(options.octaves) !== options.octaves
      ) {
        throw new Error('options.octaves must be an integer')
      }
      this.octaves = options.octaves
    } else this.octaves = 1

    if (options.hasOwnProperty('persistence')) {
      if (typeof options.persistence !== 'number') throw new Error('options.persistence must be a number')
      this.persistence = options.persistence
    } else this.persistence = 0.5

    if (options.hasOwnProperty('random')) {
      if (typeof options.random !== 'function') throw new Error('options.random must be a function')
      this.random = options.random
    } else this.random = Math.random

    let min: number
    if (options.hasOwnProperty('min')) {
      if (typeof options.min !== 'number') throw new Error('options.min must be a number')
      min = options.min
    } else min = -1

    let max: number
    if (options.hasOwnProperty('max')) {
      if (typeof options.max !== 'number') throw new Error('options.max must be a number')
      max = options.max
    } else max = 1

    if (min >= max) throw new Error(`options.min (${min}) must be less than options.max (${max})`)

    this.scale = min === -1 && max === 1
      ? value => value
      : value => min + ((value + 1) / 2) * (max - min)

    const p = new Uint8Array(256)
    for (let i = 0; i < 256; i++) p[i] = i

    let n: number
    let q: number
    for (let i = 255; i > 0; i--) {
      n = Math.floor((i + 1) * this.random())
      q = p[i]
      p[i] = p[n]
      p[n] = q
    }

    this.perm = new Uint8Array(512)
    this.permMod12 = new Uint8Array(512)
    for (let i = 0; i < 512; i++) {
      this.perm[i] = p[i & 255]
      this.permMod12[i] = this.perm[i] % 12
    }
  }

  raw2D (x: number, y: number): number {
    // Skew the input space to determine which simplex cell we're in
    const s = (x + y) * 0.3660254037844386 // Hairy factor for 2D
    const i = Math.floor(x + s)
    const j = Math.floor(y + s)
    const t = (i + j) * FastSimplexNoise.G2
    const X0 = i - t // Unskew the cell origin back to (x,y) space
    const Y0 = j - t
    const x0 = x - X0 // The x,y distances from the cell origin
    const y0 = y - Y0

    // Determine which simplex we are in.
    const i1 = x0 > y0 ? 1 : 0
    const j1 = x0 > y0 ? 0 : 1

    // Offsets for corners
    const x1 = x0 - i1 + FastSimplexNoise.G2
    const y1 = y0 - j1 + FastSimplexNoise.G2
    const x2 = x0 - 1.0 + 2.0 * FastSimplexNoise.G2
    const y2 = y0 - 1.0 + 2.0 * FastSimplexNoise.G2

    // Work out the hashed gradient indices of the three simplex corners
    const ii = i & 255
    const jj = j & 255
    const gi0 = this.permMod12[ii + this.perm[jj]]
    const gi1 = this.permMod12[ii + i1 + this.perm[jj + j1]]
    const gi2 = this.permMod12[ii + 1 + this.perm[jj + 1]]

    // Calculate the contribution from the three corners
    const t0 = 0.5 - x0 * x0 - y0 * y0
    const grad0 = FastSimplexNoise.GRAD3D[gi0]
    const n0 = t0 < 0 ? 0.0 : t0 * t0 * t0 * t0 * (grad0[0] * x0 + grad0[1] * y0)
    const t1 = 0.5 - x1 * x1 - y1 * y1
    const grad1 = FastSimplexNoise.GRAD3D[gi0]
    const n1 = t1 < 0 ? 0.0 : t1 * t1 * t1 * t1 * (grad1[0] * x1 + grad1[1] * y1)
    const t2 = 0.5 - x2 * x2 - y2 * y2
    const grad2 = FastSimplexNoise.GRAD3D[gi0]
    const n2 = t2 < 0 ? 0.0 : t2 * t2 * t2 * t2 * (grad2[0] * x2 + grad2[1] * y2)

    // Add contributions from each corner to get the final noise value.
    // The result is scaled to return values in the interval [-1, 1]
    return 70.14805770653952 * (n0 + n1 + n2)
  }

  raw3D (x: number, y: number, z: number): number {
    // Skew the input space to determine which simplex cell we're in
    const s = (x + y + z) / 3.0 // Very nice and simple skew factor for 3D
    const i = Math.floor(x + s)
    const j = Math.floor(y + s)
    const k = Math.floor(z + s)
    const t = (i + j + k) * FastSimplexNoise.G3
    const X0 = i - t // Unskew the cell origin back to (x,y,z) space
    const Y0 = j - t
    const Z0 = k - t
    const x0 = x - X0 // The x,y,z distances from the cell origin
    const y0 = y - Y0
    const z0 = z - Z0

    // Deterine which simplex we are in
    let i1: number, j1: number, k1: number // Offsets for second corner of simplex in (i,j,k) coords
    let i2: number, j2: number, k2: number // Offsets for third corner of simplex in (i,j,k) coords
    if (x0 >= y0) {
      if (y0 >= z0) {
        i1 = i2 = j2 = 1
        j1 = k1 = k2 = 0
      } else if (x0 >= z0) {
        i1 = i2 = k2 = 1
        j1 = k1 = j2 = 0
      } else {
        k1 = i2 = k2 = 1
        i1 = j1 = j2 = 0
      }
    } else {
      if (y0 < z0) {
        k1 = j2 = k2 = 1
        i1 = j1 = i2 = 0
      } else if (x0 < z0) {
        j1 = j2 = k2 = 1
        i1 = k1 = i2 = 0
      } else {
        j1 = i2 = j2 = 1
        i1 = k1 = k2 = 0
      }
    }

    const x1 = x0 - i1 + FastSimplexNoise.G3 // Offsets for second corner in (x,y,z) coords
    const y1 = y0 - j1 + FastSimplexNoise.G3
    const z1 = z0 - k1 + FastSimplexNoise.G3
    const x2 = x0 - i2 + 2.0 * FastSimplexNoise.G3 // Offsets for third corner in (x,y,z) coords
    const y2 = y0 - j2 + 2.0 * FastSimplexNoise.G3
    const z2 = z0 - k2 + 2.0 * FastSimplexNoise.G3
    const x3 = x0 - 1.0 + 3.0 * FastSimplexNoise.G3 // Offsets for last corner in (x,y,z) coords
    const y3 = y0 - 1.0 + 3.0 * FastSimplexNoise.G3
    const z3 = z0 - 1.0 + 3.0 * FastSimplexNoise.G3

    // Work out the hashed gradient indices of the four simplex corners
    const ii = i & 255
    const jj = j & 255
    const kk = k & 255
    const gi0 = this.permMod12[ii + this.perm[jj + this.perm[kk]]]
    const gi1 = this.permMod12[ii + i1 + this.perm[jj + j1 + this.perm[kk + k1]]]
    const gi2 = this.permMod12[ii + i2 + this.perm[jj + j2 + this.perm[kk + k2]]]
    const gi3 = this.permMod12[ii + 1 + this.perm[jj + 1 + this.perm[kk + 1]]]

    // Calculate the contribution from the four corners
    const t0 = 0.5 - x0 * x0 - y0 * y0 - z0 * z0
    const grad0 = FastSimplexNoise.GRAD3D[gi0]
    const n0 = t0 < 0 ? 0.0 : t0 * t0 * t0 * t0 * (grad0[0] * x0 +  grad0[1] * y0 + grad0[2] * z0)
    const t1 = 0.5 - x1 * x1 - y1 * y1 - z1 * z1
    const grad1 = FastSimplexNoise.GRAD3D[gi1]
    const n1 = t1 < 0 ? 0.0 : t1 * t1 * t1 * t1  * (grad1[0] * x1 +  grad1[1] * y1 + grad1[2] * z1)
    const t2 = 0.5 - x2 * x2 - y2 * y2 - z2 * z2
    const grad2 = FastSimplexNoise.GRAD3D[gi2]
    const n2 = t2 < 0 ? 0.0 : t2 * t2 * t2 * t2  * (grad2[0] * x2 +  grad2[1] * y2 + grad2[2] * z2)
    const t3 = 0.5 - x3 * x3 - y3 * y3 - z3 * z3
    const grad3 = FastSimplexNoise.GRAD3D[gi3]
    const n3 = t3 < 0 ? 0.0 : t3 * t3 * t3 * t3  * (grad3[0] * x3 +  grad3[1] * y3 + grad3[2] * z3)

    // Add contributions from each corner to get the final noise value.
    // The result is scaled to stay just inside [-1,1]
    return 94.68493150681972 * (n0 + n1 + n2 + n3)
  }

  scaled2D (x: number, y: number): number {
    let amplitude = this.amplitude
    let frequency = this.frequency
    let maxAmplitude = 0
    let noise = 0

    for (let i = 0; i < this.octaves; i++) {
      noise += this.raw2D(x * frequency, y * frequency) * amplitude
      maxAmplitude += amplitude
      amplitude *= this.persistence
      frequency *= 2
    }

    return this.scale(noise / maxAmplitude)
  }

  scaled3D (x: number, y: number, z: number): number {
    let amplitude = this.amplitude
    let frequency = this.frequency
    let maxAmplitude = 0
    let noise = 0

    for (let i = 0; i < this.octaves; i++) {
      noise += this.raw3D(x * frequency, y * frequency, z * frequency) * amplitude
      maxAmplitude += amplitude
      amplitude *= this.persistence
      frequency *= 2
    }

    return this.scale(noise / maxAmplitude)
  }
}

export default FastSimplexNoise;
