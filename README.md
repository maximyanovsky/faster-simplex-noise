## About this fork
Some calculations were inlined (dot function, for example) to gain performance boost when using noise function actively.
Also, spherical, cylindrical and 4D functions were removed in order to reduce code size.

# Fast Simplex Noise

[![build](https://img.shields.io/travis/joshforisha/fast-simplex-noise-js.svg)](https://travis-ci.org/joshforisha/fast-simplex-noise)
[![npm](https://img.shields.io/npm/v/fast-simplex-noise.svg)](https://www.npmjs.com/package/fast-simplex-noise)

**Important note:** _(2017-02-16)_ I have stopped supporting this library, and am instead concentrating on an OpenSimplex noise alternative, [open-simplex-noise](https://www.npmjs.com/package/open-simplex-noise). Also, see [fractal-noise](https://www.npmjs.com/package/fractal-noise) if you need fractal (octaved/shaped noise) functions.

A JavaScript implementation of the improved, faster Simplex algorithm outlined in Stefan Gustavson's [Simplex noise demystified](http://webstaff.itn.liu.se/~stegu/simplexnoise/simplexnoise.pdf). Convenience functions have been added as well, in order to avoid needing to scale and handle the raw noise data directly.

## Install

    npm install fast-simplex-noise

## Example

```js
// Generate 2D noise in a 1024x768 grid, scaled to [0, 255]

const FastSimplexNoise = require('fast-simplex-noise').default
const noiseGen = new FastSimplexNoise({ frequency: 0.01, max: 255, min: 0, octaves: 8 })

for (let x = 0; x < 1024; x++) for (let y = 0; y < 768; y++) {
  grid[x][y] = noiseGen.scaled([x, y])
}
```

### Seeded Values

You can pass a random number generator as an option (see constructor options below); specifically [seedrandom](https://www.npmjs.com/package/seedrandom) is recommended when seeded values are desired.

```js
const seedrandom = require('seedrandom')

const rng = seedrandom('hello')
const noiseGen = new FastSimplexNoise({ random: rng })
```

## API

### Constructor

#### FastSimplexNoise(options?: Options = {})

Options contains:

- `amplitude: number` – The base amplitude (default: `1.0`)
- `frequency: number` – The base frequency (default: `1.0`)
- `max: number` – The maximum scaled value to return (effective default: `1.0`)
- `min: number` – The minimum scaled value to return (effective default: `-1.0`)
- `octaves: number` – Integer; the number of octaves to sum for noise generation (default: `1`)
- `persistence: number` – The persistence of amplitude per octave (default: `0.5`)
- `random: () => number` – A function that generates random values between 0 and 1 (default: `Math.random`)

### Instance Methods

#### `raw2D(x: number, y: number): number`

Specific `raw()` call for a 2D point at (`x`, `y`).

#### `raw3D(x: number, y: number, z: number): number`

Specific `raw()` call for a 3D point at (`x`, `y`, `z`).

#### `scaled2D(x: number, y: number): number`

Specific `scaled()` call for a 2D point at (`x`, `y`).

#### `scaled3D(x: number, y: number, z: number): number`

Specific `scaled()` call for a 3D point at (`x`, `y`, `z`).

