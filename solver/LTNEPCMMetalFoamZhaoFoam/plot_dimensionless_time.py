import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent

parser = argparse.ArgumentParser(
    description="Compare runtime melt fraction with the digitized Fig. 4 benchmark."
)
parser.add_argument(
    "--metal",
    choices=("Copper", "Aluminum", "Nickel"),
    default="Copper",
    help="Metal-foam benchmark curve to overlay.",
)
parser.add_argument(
    "--data",
    type=Path,
    default=Path("postProcessing/dimensionlessTime/dimensionlessTime.dat"),
    help="Runtime dimensionless-time data file.",
)
parser.add_argument(
    "--benchmark",
    type=Path,
    default=SCRIPT_DIR / "benchmark_fig4_digitized.csv",
    help="Digitized benchmark CSV.",
)
args = parser.parse_args()

data = []
for line in args.data.read_text().splitlines():
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    values = [float(value) for value in line.split()]
    if len(values) < 9:
        raise RuntimeError(
            f"Expected 9 columns in {args.data}, found {len(values)}"
        )
    data.append(values)

if not data:
    raise RuntimeError(f"No data found in {args.data}")

benchmark = []
with args.benchmark.open(newline="") as stream:
    for row in csv.DictReader(stream):
        if row["case"] == args.metal:
            benchmark.append((float(row["tau"]), float(row["melt_fraction"])))

if not benchmark:
    raise RuntimeError(f"No {args.metal} benchmark data found in {args.benchmark}")

columns = list(zip(*data))
melt_fraction = columns[5]
tau = columns[8]

for index, value in enumerate(melt_fraction):
    if value >= 0.999:
        melt_fraction = melt_fraction[: index + 1]
        tau = tau[: index + 1]
        break

colors = {"Copper": "black", "Aluminum": "blue", "Nickel": "red"}
color = colors[args.metal]

plt.plot(
    [point[0] for point in benchmark],
    [point[1] for point in benchmark],
    color=color,
    linestyle="--",
    linewidth=2,
    label=f"{args.metal} paper benchmark",
)
plt.plot(
    tau,
    melt_fraction,
    color=color,
    linewidth=3,
    label=f"{args.metal} numerical",
)
plt.xlabel(r"$\tau = Ste\,Fo$")
plt.ylabel("Melt fraction")
plt.xlim(left=0)
plt.ylim(0, 1.02)
plt.grid(alpha=0.25)
plt.legend()
plt.tight_layout()

output = Path(f"melt_fraction_{args.metal.lower()}_numerical_vs_benchmark.png")
plt.savefig(output, dpi=300)
print(f"Wrote {output}")
