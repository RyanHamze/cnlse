from cnlse.app import build_simulation


def test_build_simulation_returns_matching_traces():
    data = build_simulation(
        symbols=8,
        samples_per_symbol=24,
        symbol_rate_gbaud=10,
        pattern_type="De Bruijn",
        seed=7,
        phi_rad=0.5,
        extinction_ratio_db=13,
        wavelength_nm=1550,
        tx_length_km=10,
        tx_loss_db_km=0.2,
        tx_dispersion=17,
        comp_loss_db_km=0.6,
        comp_dispersion=-100,
        nonlinear=False,
    )

    expected_samples = 8 * 24
    assert len(data["bits"]) == 8
    assert data["time"].shape == (expected_samples,)
    assert data["input_power"].shape == (expected_samples,)
    assert data["uncompensated_power"].shape == (expected_samples,)
    assert data["compensated_power"].shape == (expected_samples,)
    assert data["comp_length_km"] > 0
