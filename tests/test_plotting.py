import numpy as np

from cnlse.plotting import plot_power_traces


def test_plot_power_traces_writes_file(tmp_path):
    output = tmp_path / "plot.png"
    time = np.linspace(0, 1, 16)
    result = plot_power_traces(time, [("demo", time**2)], output)
    assert result == output
    assert output.exists()
    assert output.stat().st_size > 0
