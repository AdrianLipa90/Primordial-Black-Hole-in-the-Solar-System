#!/usr/bin/env python3
from __future__ import annotations
import math

RA_DEG = 107.576646
DEC_DEG = -25.907579
GAL_L_DEG = 238.1709218507047
GAL_B_DEG = -7.665782897020605


def _spread_bits(value: int) -> int:
    out = 0
    bit = 0
    while value:
        if value & 1:
            out |= 1 << (2 * bit)
        value >>= 1
        bit += 1
    return out


def ang2pix_nest(nside: int, lon_deg: float, lat_deg: float) -> int:
    if nside <= 0 or nside & (nside - 1):
        raise ValueError("nside must be a positive power of two")
    phi = math.radians(lon_deg) % (2.0 * math.pi)
    z = math.sin(math.radians(lat_deg))
    za = abs(z)
    tt = phi / (0.5 * math.pi)

    if za <= 2.0 / 3.0:
        temp1 = nside * (0.5 + tt)
        temp2 = nside * (0.75 * z)
        jp = math.floor(temp1 - temp2)
        jm = math.floor(temp1 + temp2)
        ifp = jp // nside
        ifm = jm // nside
        if ifp == ifm:
            face = ifp | 4
        elif ifp < ifm:
            face = ifp
        else:
            face = ifm + 8
        ix = jm % nside
        iy = nside - (jp % nside) - 1
    else:
        ntt = math.floor(tt)
        tp = tt - ntt
        tmp = nside * math.sqrt(3.0 * (1.0 - za))
        jp = min(math.floor(tp * tmp), nside - 1)
        jm = min(math.floor((1.0 - tp) * tmp), nside - 1)
        if z >= 0:
            face = ntt
            ix = nside - jm - 1
            iy = nside - jp - 1
        else:
            face = ntt + 8
            ix = jp
            iy = jm

    ipf = _spread_bits(ix) | (_spread_bits(iy) << 1)
    return int(face * nside * nside + ipf)


def main() -> None:
    p8 = ang2pix_nest(8, GAL_L_DEG, GAL_B_DEG)
    p512 = ang2pix_nest(512, GAL_L_DEG, GAL_B_DEG)
    p2048 = ang2pix_nest(2048, GAL_L_DEG, GAL_B_DEG)

    assert p8 == 488
    assert p512 >> (2 * (9 - 3)) == p8
    assert p2048 >> (2 * (11 - 3)) == p8

    print({
        "order3_tile": p8,
        "nside512_pixel": p512,
        "nside2048_pixel": p2048,
        "relative_tile_path": f"Norder3/Dir{(p8 // 10000) * 10000}/Npix{p8}.fits",
    })


if __name__ == "__main__":
    main()
