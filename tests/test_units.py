import numpy as np

from v1lb.units import amplitude_cutoff, isi_contamination, presence_ratio


def test_isi_contamination_zero_for_refractory_train_and_grows_with_mixing():
    rng = np.random.default_rng(0)
    clean = np.cumsum(0.002 + rng.exponential(0.1, 5000))          # 2 ms dead time: no violations
    assert isi_contamination(clean, clean[-1]) == 0
    other = rng.uniform(0, clean[-1], 2500)                         # a third of spikes from a Poisson neighbour
    mixed = np.sort(np.concatenate([clean, other]))
    assert 0.1 < isi_contamination(mixed, clean[-1]) < 1.0


def test_amplitude_cutoff_detects_truncation():
    rng = np.random.default_rng(1)
    amp = rng.normal(20, 4, 20000)
    assert amplitude_cutoff(amp) < 0.02
    assert 0.3 < amplitude_cutoff(amp[amp > 20]) <= 0.5                # half the distribution below threshold


def test_presence_ratio():
    s = np.concatenate([np.arange(0, 300, 0.5), np.arange(540, 600, 0.5)])   # silent 300-540 s
    assert np.isclose(presence_ratio(s, 0, 600), 6 / 10)


def test_base_ticks_maps_joined_packets():
    from types import SimpleNamespace
    from v1lb.units import PhySort, base_ticks
    ns6 = SimpleNamespace(t_start=102, period=1)
    s = PhySort(None, None, np.array([0, 999, 1000, 1500]), np.zeros(4, int), np.ones(4), None,
                packets=[[102, 1000], [1130, 5000]])                   # 28-tick gap between packets
    assert base_ticks(s, ns6).tolist() == [102, 1101, 1130, 1630]
    s.packets = None
    assert base_ticks(s, ns6).tolist() == [102, 1101, 1102, 1602]
