import numpy as np

from v1lb.binocular import CONTRASTS, equivalent_contrast, eye_weight, invert_crf


def nr(c, rmax=10.0, c50=0.15, n=2.5):
    c = np.asarray(c, float)
    return rmax * c ** n / (c ** n + c50 ** n)


def test_invert_crf_recovers_contrast():
    r = nr(CONTRASTS)
    c, cens = invert_crf(np.repeat(r[:, None], 3, axis=1), np.array([nr(0.08), nr(0.3), nr(0.7)]))
    assert np.allclose(c, [0.08, 0.3, 0.7], rtol=0.15) and (cens == 0).all()
    _, cens = invert_crf(r, np.array(nr(0.9) * 1.2))
    assert cens == 1


def test_eye_weight_recovers_scale():
    a = 0.6
    assert abs(eye_weight(nr(CONTRASTS), nr(a * CONTRASTS)) - a) < 0.08


def test_k_is_latency_proof_and_separates_rules():
    t = np.arange(0, 200)

    def resp(c_eff):                       # higher effective contrast: larger and earlier (latency 60 -> 40 ms)
        lat = 40 + 20 * (1 - np.clip(np.log10(np.maximum(c_eff, 1e-3) / 0.055) / np.log10(0.9 / 0.055), 0, 1))
        return nr(c_eff) * np.clip((t - lat) / 10, 0, 1)

    mono = np.stack([resp(c) for c in CONTRASTS])                       # [contrast, time]
    a, c = 1.0, 0.11
    lin = resp((1 + a) * c)                                              # linear summation, earlier onset
    inv = resp(c)                                                        # ocularity invariance
    on = t >= 60
    k_lin, _ = equivalent_contrast(mono, lin, c)
    k_inv, _ = equivalent_contrast(mono, inv, c)
    assert np.allclose(k_lin[on], 2.0, atol=0.1) and np.allclose(k_inv[on], 1.0, atol=0.05)
    # a pure latency difference (binocular = monocular at 2c, including its earlier onset) keeps k constant
    assert np.nanstd(k_lin[on]) < 0.05
