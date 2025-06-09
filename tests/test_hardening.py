import math

# Reference implementations matching docs/CMLformat07.md

def H_iso(kappa, yld, hk, hpa, hpb):
    return yld + hk * kappa + (hpa - yld) * (1 - math.exp(-hpb * kappa))

def dH_iso_dk(kappa, yld, hk, hpa, hpb):
    return hk + (hpa - yld) * hpb * math.exp(-hpb * kappa)

def H_kin(kappa, hk, theta):
    return (1 - theta) * hk * kappa

def dH_kin_dk(hk, theta):
    return (1 - theta) * hk

def test_values():
    kappa = 0.2
    yld, hk, hpa, hpb, theta = 100.0, 50.0, 120.0, 10.0, 0.3
    iso = H_iso(kappa, yld, hk, hpa, hpb)
    kin = H_kin(kappa, hk, theta)
    assert abs(iso - (yld + hk*kappa + (hpa - yld)*(1-math.exp(-hpb*kappa)))) < 1e-12
    assert abs(dH_iso_dk(kappa, yld, hk, hpa, hpb) - (hk + (hpa - yld)*hpb*math.exp(-hpb*kappa))) < 1e-12
    assert abs(kin - ((1-theta)*hk*kappa)) < 1e-12
    assert abs(dH_kin_dk(hk, theta) - ((1-theta)*hk)) < 1e-12

if __name__ == '__main__':
    test_values()
    print('hardening tests passed')
