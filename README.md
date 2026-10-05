# Spacecraft Telecom Link Budget Calculator

A small Python toolkit for estimating the RF link budget between a spacecraft and an Earth ground station. Given the spacecraft, ground station, orbit and payload parameters, it computes antenna gains, EIRP, path loss, atmospheric attenuation (via the ITU-R models), G/T, C/N0, Eb/N0 and the final **link margin**, for both uplink and downlink.

## Project structure

```
.
├── payload.py   # Payload class: imaging payload parameters used to derive the data rate
└── telecom.py   # Downlink class: link budget calculations
```

| File | Purpose |
|------|---------|
| `payload.py` | Holds the payload description (swath width, pixel size, bits per pixel, duty cycle, downlink time). |
| `telecom.py` | Contains the `Downlink` class, which performs all link budget calculations and takes a `Payload` instance. |

## Installation

### 1. Prerequisites

- Python 3.8 or newer
- `pip` (or conda)

### 2. (Recommended) Create a virtual environment

```bash
python -m venv venv
source venv/bin/activate        # Linux / macOS
venv\Scripts\activate           # Windows
```

### 3. Install the dependencies

`telecom.py` needs `numpy` and [`itur`](https://github.com/inigodelportillo/ITU-Rpy), a Python implementation of the ITU-R propagation recommendations. It is used here for gaseous attenuation on the slant path (ITU-R P.676).

```bash
pip install numpy itur
```

`itur` pulls in its own dependencies automatically (e.g. `scipy`, `astropy`, `pyproj`, `matplotlib`).

**Using conda instead:**

```bash
conda install -c conda-forge itur
```

**Verify the installation:**

```bash
python -c "import itur; print(itur.__version__)"
```

If this prints a version number, you are ready to go.

## Usage

Place `payload.py` and `telecom.py` in the same folder (`telecom.py` does `from payload import Payload`).

```python
from payload import Payload
from telecom import Downlink

# 1. Describe the payload
payload = Payload(
    swath_width_angle=10,   # deg
    pixel_size=0.1,         # arcmin
    bits_per_pixel=8,
    duty_cycle=0.3,         # fraction of the orbit the instrument is on
    downlink_time=4,        # hours per day available for downlink
)

# 2. Describe the link
link = Downlink(
    f_downlink=8.2,              # GHz
    sc_tx_power=10,              # W
    gs_tx_power=100,             # W
    sc_antenna_diameter=0.5,     # m
    gs_antenna_diameter=3.0,     # m
    turn_around_ratio=221 / 240,
    orbit_altitude=500,          # km
    elongation_angle=0,          # deg (only used for interplanetary missions)
    pointing_offset_angle=0.1,   # deg
    required_uplink_dr=2000,     # bps
    payload=payload,
    required_bit_error_rate=1e-5,
    elevation_angle=10,          # deg (default)
)

# 3. Compute the link margin
MU_EARTH = 3.986004418e14        # m^3/s^2
R_EARTH = 6371.0                 # km

margin, status = link.calculate_link_margin(
    mode="downlink",             # or "uplink"
    type="Earth Orbit",          # "Earth Orbit", "Lunar Orbit", or anything else for interplanetary
    d_S=149597870.7,             # km, spacecraft/planet to Sun distance (interplanetary only)
    mu=MU_EARTH,
    body_radius=R_EARTH,
)

print(f"Link margin: {margin:.2f} dB")
print(status)
```

### Mission types

The `type` argument selects how the slant range is computed:

| `type` | Slant range used |
|--------|------------------|
| `"Earth Orbit"` | Geometric slant range from altitude and minimum elevation angle |
| `"Lunar Orbit"` | Earth-Moon distance (384,400 km) |
| anything else | Interplanetary: law of cosines using the Earth-Sun distance (1 AU), `d_S` and the elongation angle |

## Main methods (`Downlink`)

| Method | Returns |
|--------|---------|
| `determine_pointing_loss(mode)` | Pointing loss in dB |
| `determine_free_space_path_loss(mode, type, d_S)` | Free space path loss in dB |
| `determine_atmospheric_attenuation(mode)` | Gaseous attenuation in dB (ITU-R, clear sky) |
| `calculate_transmission_gain(mode)` / `calculate_reception_gain(mode)` | Antenna gains in dBi |
| `calculate_eirp(mode)` | EIRP in dBW |
| `calculate_g_over_t(mode, type)` | Receiver G/T in dB/K |
| `calculate_c_over_n0(mode, type, d_S)` | Carrier-to-noise density in dB-Hz |
| `determine_data_rate(mu, body_radius, mode)` | Required data rate in bps |
| `calculate_eb_over_n0(mode, type, d_S, body_radius, mu)` | Eb/N0 in dB |
| `calculate_link_margin(mode, type, d_S, mu, body_radius, eb_n0_threshold_dB=3.0)` | Link margin in dB and a pass/fail message |

`mode` is `"downlink"` or `"uplink"`. The uplink frequency is `f_downlink * turn_around_ratio`.

## Units

| Parameter | Unit |
|-----------|------|
| Frequencies | GHz |
| Transmit power | W |
| Antenna diameters | m |
| Orbit altitude, body radius, `d_S` | km |
| `mu` (gravitational parameter) | m³/s² |
| Angles (elevation, elongation, pointing offset, swath width) | degrees |
| Payload `pixel_size` | arcminutes |
| Payload `downlink_time` | hours per day |
| Data rates | bits per second |

## Assumptions and limitations

- Both antennas are parabolic with a fixed efficiency of 0.55.
- Atmospheric attenuation uses only gaseous absorption (clear sky) at standard sea-level conditions; no rain, cloud or scintillation effects, and no ground station coordinates.
- The elevation angle used for the ITU-R model is bounded to a minimum of 5° (model validity).
- Antenna temperature on the uplink is assumed to be 290 K; on the downlink it is derived from ITU-R P.618 with `T_mr = 275 K`.
- Circular orbit is assumed for the payload data-rate calculation.
- Link margin is computed against an Eb/N0 threshold of 3 dB by default (per ECSS-E-ST-50-05C).
- `required_bit_error_rate` is stored but not currently used in the margin calculation.
- The Earth-orbit slant range ignores Earth's rotation and ground station altitude.

## References

- ITU-R P.618, P.676 (propagation data and prediction methods)
- T. Pratt and J. Allnutt, *Satellite Communications*
- ECSS-E-ST-50-05C, Radio frequency and modulation
- [ITU-Rpy documentation](https://itur.readthedocs.io/)
