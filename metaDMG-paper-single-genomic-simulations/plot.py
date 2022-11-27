import matplotlib.pyplot as plt
import numpy as np

plt.style.use("plotstyle.txt")

#%%


t = np.linspace(0.0, 1.0, 100)


fig, ax = plt.subplots(
    # figsize=(6, 4),
    # tight_layout=True,
)

for i in range(10):
    s = np.cos(4 * np.pi * t) + i
    ax.plot(t, s, label=i)

ax.set_xlabel(
    # r"\textbf{time (s)}",
    r"time (s)",
)
ax.set_ylabel(
    "\\textit{Velocity (\N{DEGREE SIGN}/sec)}",
    # fontsize=16,
)
ax.set_title(
    r"ABC123 vs $\mathrm{ABC123}^{123}$"
    # "Here a long text title",
    # r"\TeX\ is Number $\displaystyle\sum_{n=1}^\infty" r"\frac{-e^{i\pi}}{2^n}$!",
    # fontsize=16,
    # color="r",
)

ax.legend()

fig.savefig("test-fig.pdf")
