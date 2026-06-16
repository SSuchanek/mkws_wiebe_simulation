"""
Symulacja termodynamiczna silnika tłokowego z modelem spalania Wiebescha
oraz analizą kinetyki chemicznej spalania w Cantera
=========================================================================
Przedmiot: Metody Komputerowe w Spalaniu
Autor: [Twoje imię i nazwisko]

Model obejmuje:
  - Analizę kinetyki chemicznej spalania CH4/powietrze (Cantera, mech. GRI-3.0):
      * adiabatyczna temperatura spalania T_ad (równowaga HP)
      * równowagowy skład spalin
      * prędkość spalania laminarnego S_L (model płomienia 1D, FreeFlame)
  - Geometrię układu korbowo-tłokowego
  - Sprężanie i rozprężanie (gaz doskonały, γ = const)
  - Spalanie w cylindrze opisane funkcją Wiebescha
  - Bilans energii (I zasada termodynamiki)
  - Wykresy: p-V, p(φ), T(φ), x(φ), T_ad(λ), S_L(λ)
  - Wskaźniki: p_max, T_max, T_ad, S_L, IMEP, η_th
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cantera as ct

# ─────────────────────────────────────────────────────────────
# 1. PARAMETRY SILNIKA
# ─────────────────────────────────────────────────────────────

# Typ silnika: 'ZI' (Otto / iskrowy) lub 'ZS' (Diesel / samoczynny)
TYP_SILNIKA = 'ZI'

# Geometria
eps   = 10.0        # [-]   stopień sprężania
Vd    = 500e-6      # [m³]  pojemność skokowa (500 cm³)
r_kl  = 0.5        # [-]   stosunek ramienia korby do długości korbowodu (r/l)

# Warunki dolotowe
p1    = 1.01325e5   # [Pa]  ciśnienie na początku sprężania
T1    = 320.0       # [K]   temperatura na początku sprężania

# Skład mieszanki
lam   = 1.0         # [-]   współczynnik nadmiaru powietrza (λ=1 → stechiometria)

# ─────────────────────────────────────────────────────────────
# 2. PARAMETRY MODELU WIEBESCHA
# ─────────────────────────────────────────────────────────────

m        = 2.0      # [-]     wykładnik kształtu (ZI≈2, ZS≈0.8)
phi0     = -10.0    # [°OWK]  kąt początku spalania (przed GMP < 0)
delta_phi = 60.0    # [°OWK]  czas trwania (kąt) spalania

# Sprawność wskazana (uwzględnia straty cieplne)
eta_i    = 0.95     # [-]

# ─────────────────────────────────────────────────────────────
# 3. STAŁE FIZYCZNE I PALIWO
# ─────────────────────────────────────────────────────────────

R_air  = 287.0      # [J/(kg·K)]  stała gazowa powietrza
gamma  = 1.38       # [-]         wykładnik adiabaty (uśredniony)
Cv     = R_air / (gamma - 1)
Cp     = gamma * Cv

if TYP_SILNIKA == 'ZI':
    Hu   = 50.0e6   # [J/kg]  wartość opałowa CH4 (zgodna z analizą Cantera, sekcja 3a)
    L_st = 17.2     # [-]     stechiometryczne zapotrzebowanie powietrza dla CH4 (masowe)
else:
    Hu   = 42.5e6   # [J/kg]  wartość opałowa oleju napędowego (model ZS — tabelaryczny)
    L_st = 14.5

# ─────────────────────────────────────────────────────────────
# 3a. KINETYKA CHEMICZNA SPALANIA — CANTERA (mechanizm GRI-3.0)
# ─────────────────────────────────────────────────────────────
#
# Paliwo modelowe: CH4 (metan) — substytut paliwa rzeczywistego,
# ponieważ mechanizm GRI-3.0 nie zawiera węglowodorów ciężkich
# (typu n-oktan/dodekan stosowanych jako surogaty benzyny/oleju
# napędowego). Metan jest standardowym paliwem referencyjnym
# do analizy kinetyki spalania w Cantera.
#
# Liczone wielkości:
#   - T_ad   : adiabatyczna temperatura spalania (równowaga H=const, p=const)
#   - skład spalin w równowadze chemicznej
#   - S_L    : prędkość spalania laminarnego (model płomienia 1D, FreeFlame)
#
# Wyniki wykorzystywane są jako fizyczna walidacja modelu Wiebescha
# (sekcja 6) — temperatura T_ad wyznacza górną granicę temperatury
# osiąganej w cylindrze.

PALIWO_MECH = 'gri30.yaml'   # mechanizm reakcji chemicznych (53 gatunki, 325 reakcji)
PALIWO_NAZWA = 'CH4'         # paliwo modelowe (metan)
OKSYDATOR = 'O2:1, N2:3.76'  # powietrze (uproszczone, bez Ar)


def analiza_spalania_cantera(lam, T1, p1, mech=PALIWO_MECH, paliwo=PALIWO_NAZWA):
    """
    Analiza spalania w Cantera: temperatura adiabatyczna, skład
    spalin w równowadze oraz prędkość spalania laminarnego.

    Parametry
    ---------
    lam : współczynnik nadmiaru powietrza [-]  (λ = 1/φ_eq)
    T1  : temperatura mieszanki przed spalaniem [K]
    p1  : ciśnienie mieszanki przed spalaniem [Pa]

    Zwraca
    ------
    słownik: T_ad [K], skład_spalin (dict mol. udziałów > 0.1%), S_L [m/s]
    """
    phi_eq = 1.0 / lam   # współczynnik ekwiwalencji (Cantera convention)

    # ── Adiabatyczna temperatura spalania i równowaga chemiczna ──
    gas = ct.Solution(mech)
    gas.set_equivalence_ratio(phi_eq, paliwo, OKSYDATOR)
    gas.TP = T1, p1
    gas.equilibrate('HP')   # stałe H i p -> adiabatyczne spalanie izobaryczne

    T_ad = gas.T
    sklad_spalin = {
        sp: gas[sp].X[0]
        for sp in gas.species_names
        if gas[sp].X[0] > 1e-3
    }

    # ── Prędkość spalania laminarnego (model płomienia 1D) ──
    gas_flame = ct.Solution(mech)
    gas_flame.set_equivalence_ratio(phi_eq, paliwo, OKSYDATOR)
    gas_flame.TP = T1, p1
    flame = ct.FreeFlame(gas_flame, width=0.03)
    flame.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
    try:
        flame.solve(loglevel=0, auto=True)
        S_L = flame.velocity[0]   # [m/s]
    except ct.CanteraError:
        S_L = np.nan   # poza zakresem zapłonu (np. bardzo bogata/biedna mieszanka)

    return {'T_ad': T_ad, 'sklad_spalin': sklad_spalin, 'S_L': S_L}


# Analiza spalania dla aktualnego λ (parametr z sekcji 1)
wynik_cantera = analiza_spalania_cantera(lam=lam, T1=T1, p1=p1)

# ─────────────────────────────────────────────────────────────
# 4. GEOMETRIA UKŁADU KORBOWO-TŁOKOWEGO
# ─────────────────────────────────────────────────────────────

Vc  = Vd / (eps - 1)        # [m³] objętość komory spalania
V1  = Vd + Vc               # [m³] objętość max (BMP)
l_kl = 1.0 / r_kl          # [-]  unormowana długość korbowodu (r=1)

def objetosc(phi_deg: np.ndarray) -> np.ndarray:
    """
    Objętość cylindra w funkcji kąta obrotu wału korbowego.

    Parametry
    ---------
    phi_deg : kąt OWK w stopniach (0° = GMP — górne martwe położenie)

    Zwraca
    ------
    V : objętość [m³]
    """
    phi = np.deg2rad(phi_deg)
    # Przemieszczenie tłoka względem GMP (unormowane do r=1):
    # s = r*(1 - cos(φ)) + l*(1 - sqrt(1 - (r/l·sin(φ))²))
    delta = np.sqrt(l_kl**2 - np.sin(phi)**2)
    s_norm = (1 - np.cos(phi)) + l_kl * (1 - delta)
    # Rzeczywista objętość: V = Vc + (Vd/2) * (1 - cos) ... uproszczenie
    # Pełna kinematyka:
    return Vc + (Vd / 2) * s_norm / (1 + l_kl - 1)   # skalowanie do Vd

def objetosc_dokladna(phi_deg: np.ndarray) -> np.ndarray:
    """Dokładna objętość z pełną kinemayką korbowodu."""
    phi = np.deg2rad(phi_deg)
    lam_r = r_kl          # r/l
    pos = (1 - np.cos(phi)) + (1/lam_r) * (1 - np.sqrt(1 - (lam_r * np.sin(phi))**2))
    pos_max = 2.0          # przy φ=180° (BMP)
    return Vc + Vd * pos / pos_max

# ─────────────────────────────────────────────────────────────
# 5. FUNKCJA WIEBESCHA
# ─────────────────────────────────────────────────────────────

def wiebe(phi_deg: np.ndarray, phi0: float, delta_phi: float, m: float) -> np.ndarray:
    """
    Funkcja Wiebescha — udział masy spalonego paliwa.

    x(φ) = 1 - exp(-a · ((φ - φ₀) / Δφ)^(m+1))

    gdzie a = 6.908 (przyjęte η_v = 0.999)

    Parametry
    ---------
    phi_deg   : kąt OWK [°]
    phi0      : kąt początku spalania [°OWK]
    delta_phi : czas trwania spalania [°OWK]
    m         : wykładnik kształtu Wiebescha [-]

    Zwraca
    ------
    x : udział masy paliwa spalonego [0..1]
    """
    a = 6.908  # -ln(1 - 0.999)
    x = np.zeros_like(phi_deg, dtype=float)
    mask = (phi_deg > phi0) & (phi_deg < phi0 + delta_phi)
    u = (phi_deg[mask] - phi0) / delta_phi
    x[mask] = 1.0 - np.exp(-a * u**(m + 1))
    x[phi_deg >= phi0 + delta_phi] = 1.0
    return x

# ─────────────────────────────────────────────────────────────
# 6. SYMULACJA OBIEGU TERMODYNAMICZNEGO
# ─────────────────────────────────────────────────────────────

def symuluj_obieg(eps, p1, T1, lam, m, phi0, delta_phi, eta_i,
                  Hu, L_st, Vd, Vc, gamma, R_air, r_kl, dphi=0.5):
    """
    Symulacja obiegu pracy silnika tłokowego metodą różnic skończonych.

    Obieg (kąt OWK):
      -360° → -180°  dolot (pominięty — przyjmujemy p=p1, T=T1)
       -180° → 0°    sprężanie
          0° → 180°  spalanie + rozprężanie
        180° → 360°  wydech (pominięty)

    Równanie energii (I zasada dla układu zamkniętego):
      dU = δQ - δW
      m·Cv·dT = dQ_chemiczne - p·dV
      dp = (γ-1)/V · dQ - γ·p/V · dV

    Parametry
    ---------
    dphi : krok całkowania [°OWK]

    Zwraca
    ------
    słownik z tablicami wyników i wskaźnikami
    """
    phi_arr = np.arange(-180, 180 + dphi, dphi)
    V_arr   = objetosc_dokladna(phi_arr)

    # Masa czynnika roboczego
    masa_powietrza = p1 * (Vd + Vc) / (R_air * T1)
    masa_paliwa    = masa_powietrza / (lam * L_st)
    Q_total        = masa_paliwa * Hu * eta_i   # [J] ciepło doprowadzone

    # Funkcja Wiebescha — udział spalania w każdym kroku
    x_arr  = wiebe(phi_arr, phi0, delta_phi, m)
    dQ_arr = np.diff(x_arr, prepend=x_arr[0]) * Q_total  # przyrosty ciepła

    # Całkowanie
    p_arr = np.zeros(len(phi_arr))
    T_arr = np.zeros(len(phi_arr))
    p_arr[0] = p1
    T_arr[0] = T1

    masa_total = masa_powietrza + masa_paliwa

    for i in range(1, len(phi_arr)):
        dV = V_arr[i] - V_arr[i - 1]
        dQ = dQ_arr[i]
        V_i = V_arr[i]
        p_prev = p_arr[i - 1]

        # Zmiana ciśnienia (I zasada — gaz doskonały):
        # dp = (γ-1)/V · dQ - γ·p/V · dV
        dp = ((gamma - 1) / V_i) * dQ - (gamma * p_prev / V_i) * dV
        p_arr[i] = p_prev + dp
        T_arr[i] = p_arr[i] * V_arr[i] / (masa_total * R_air)

    # ── Wskaźniki ──────────────────────────────────────────────
    p_max = np.max(p_arr)
    T_max = np.max(T_arr)

    # IMEP — średnie ciśnienie indykowane
    # W = ∮ p dV  (trapezowa reguła)
    praca = np.trapezoid(p_arr, V_arr)   # [J]
    IMEP  = praca / Vd               # [Pa]

    # Sprawność termodynamiczna
    eta_th = praca / Q_total if Q_total > 0 else 0.0

    return {
        'phi': phi_arr,
        'p':   p_arr,
        'T':   T_arr,
        'V':   V_arr,
        'x':   x_arr,
        'dQ':  dQ_arr,
        'p_max': p_max,
        'T_max': T_max,
        'IMEP':  IMEP,
        'eta_th': eta_th,
        'Q_total': Q_total,
        'praca': praca,
    }

# ─────────────────────────────────────────────────────────────
# 7. URUCHOMIENIE SYMULACJI
# ─────────────────────────────────────────────────────────────

wyniki = symuluj_obieg(
    eps=eps, p1=p1, T1=T1, lam=lam,
    m=m, phi0=phi0, delta_phi=delta_phi, eta_i=eta_i,
    Hu=Hu, L_st=L_st, Vd=Vd, Vc=Vc,
    gamma=gamma, R_air=R_air, r_kl=r_kl,
    dphi=0.5
)

# ─────────────────────────────────────────────────────────────
# 8. WYDRUK WSKAŹNIKÓW
# ─────────────────────────────────────────────────────────────

print("=" * 55)
print(f"  SYMULACJA SILNIKA TŁOKOWEGO — MODEL WIEBESCHA")
print(f"  Typ silnika : {TYP_SILNIKA}")
print("=" * 55)
print(f"  Stopień sprężania ε       : {eps:.1f}")
print(f"  Ciśnienie dolotowe p₁     : {p1/1e5:.3f} bar")
print(f"  Temperatura dolotowa T₁   : {T1:.1f} K")
print(f"  Nadmiar powietrza λ        : {lam:.2f}")
print(f"  Wykładnik Wiebescha m     : {m:.1f}")
print(f"  Kąt początku spalania φ₀  : {phi0:.1f} °OWK")
print(f"  Czas spalania Δφ          : {delta_phi:.1f} °OWK")
print(f"  Sprawność wskazana η_i    : {eta_i*100:.0f} %")
print("-" * 55)
print(f"  --- Kinetyka chemiczna (Cantera, {PALIWO_NAZWA}/powietrze) ---")
print(f"  Adiabatyczna temp. spalania T_ad : {wynik_cantera['T_ad']:.1f} K")
print(f"  Prędkość spalania laminarnego S_L : {wynik_cantera['S_L']*100:.2f} cm/s")
print("  Skład spalin (równowaga, udział molowy > 0.1%):")
for sp, frac in sorted(wynik_cantera['sklad_spalin'].items(), key=lambda kv: -kv[1]):
    print(f"    {sp:>5} : {frac:.4f}")
print("-" * 55)
print(f"  Ciepło doprowadzone Q     : {wyniki['Q_total']:.2f} J")
print(f"  Praca indykowana W        : {wyniki['praca']:.2f} J")
print(f"  Ciśnienie max p_max       : {wyniki['p_max']/1e5:.2f} bar")
print(f"  Temperatura max T_max     : {wyniki['T_max']:.1f} K")
print(f"  IMEP                      : {wyniki['IMEP']/1e5:.3f} bar")
print(f"  Sprawność termiczna η_th  : {wyniki['eta_th']*100:.2f} %")
print("=" * 55)
print("  Uwaga: T_max w cylindrze > T_ad, ponieważ spalanie w cylindrze")
print("  zachodzi w mieszance już sprężonej (wysokie T, p), podczas gdy")
print("  T_ad z Cantera dotyczy spalania mieszanki o warunkach początkowych T1, p1.")

# ─────────────────────────────────────────────────────────────
# 9. WYKRESY
# ─────────────────────────────────────────────────────────────

phi = wyniki['phi']
p   = wyniki['p'] / 1e5          # bar
T   = wyniki['T']                 # K
V   = wyniki['V'] * 1e6          # cm³
x   = wyniki['x']

fig = plt.figure(figsize=(14, 14))
fig.suptitle(
    f"Symulacja termodynamiczna silnika {TYP_SILNIKA}  |  "
    f"ε = {eps}, λ = {lam}, m = {m}, φ₀ = {phi0}°, Δφ = {delta_phi}°",
    fontsize=13, fontweight='bold'
)

gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.32)

# ── 9.1  Wykres p-V ─────────────────────────────────────────
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(V, p, color='#1a6bbd', linewidth=2)
ax1.fill_between(V, p, alpha=0.10, color='#1a6bbd')
ax1.set_xlabel('Objętość V [cm³]', fontsize=11)
ax1.set_ylabel('Ciśnienie p [bar]', fontsize=11)
ax1.set_title('Wykres p-V (indykatorowy)', fontsize=11)
ax1.annotate(f'p_max = {wyniki["p_max"]/1e5:.1f} bar',
             xy=(V[np.argmax(p)], np.max(p)),
             xytext=(V[np.argmax(p)] + 5, np.max(p) * 0.9),
             arrowprops=dict(arrowstyle='->', color='gray'),
             fontsize=9, color='#333')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(left=0)
ax1.set_ylim(bottom=0)

# ── 9.2  Ciśnienie p(φ) ─────────────────────────────────────
ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(phi, p, color='#c94a1a', linewidth=2)
ax2.fill_between(phi, p, alpha=0.10, color='#c94a1a')
ax2.axvline(0, color='gray', linestyle='--', linewidth=0.8, label='GMP')
ax2.axvline(phi0, color='#888', linestyle=':', linewidth=1, label=f'φ₀ = {phi0}°')
ax2.set_xlabel('Kąt OWK φ [°]', fontsize=11)
ax2.set_ylabel('Ciśnienie p [bar]', fontsize=11)
ax2.set_title('Ciśnienie w cylindrze p(φ)', fontsize=11)
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)
ax2.set_ylim(bottom=0)

# ── 9.3  Temperatura T(φ) ───────────────────────────────────
ax3 = fig.add_subplot(gs[1, 0])
ax3.plot(phi, T, color='#1d9e75', linewidth=2)
ax3.fill_between(phi, T, alpha=0.10, color='#1d9e75')
ax3.axvline(0, color='gray', linestyle='--', linewidth=0.8, label='GMP')
ax3.axhline(T1, color='#aaa', linestyle=':', linewidth=1, label=f'T₁ = {T1} K')
ax3.set_xlabel('Kąt OWK φ [°]', fontsize=11)
ax3.set_ylabel('Temperatura T [K]', fontsize=11)
ax3.set_title('Temperatura gazów T(φ)', fontsize=11)
ax3.annotate(f'T_max = {wyniki["T_max"]:.0f} K',
             xy=(phi[np.argmax(T)], np.max(T)),
             xytext=(phi[np.argmax(T)] + 10, np.max(T) * 0.95),
             arrowprops=dict(arrowstyle='->', color='gray'),
             fontsize=9, color='#333')
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3)

# ── 9.4  Funkcja Wiebescha x(φ) ────────────────────────────
ax4 = fig.add_subplot(gs[1, 1])
phi_w = np.linspace(phi0 - 15, phi0 + delta_phi + 15, 500)
x_w   = wiebe(phi_w, phi0, delta_phi, m)
ax4.plot(phi_w, x_w, color='#7f52b8', linewidth=2.5, label=f'm = {m}')

# Pokaż też m=0.5 i m=4 dla porównania
for m_cmp, ls, alpha in [(0.5, '--', 0.5), (4.0, ':', 0.5)]:
    ax4.plot(phi_w, wiebe(phi_w, phi0, delta_phi, m_cmp),
             color='#7f52b8', linewidth=1.5, linestyle=ls, alpha=alpha,
             label=f'm = {m_cmp}')

ax4.axhline(0.5, color='gray', linestyle=':', linewidth=0.8, label='x = 0.5')
ax4.set_xlabel('Kąt OWK φ [°]', fontsize=11)
ax4.set_ylabel('Udział masy spalonej x [-]', fontsize=11)
ax4.set_title('Funkcja Wiebescha x(φ)', fontsize=11)
ax4.legend(fontsize=9)
ax4.grid(True, alpha=0.3)
ax4.set_ylim(-0.05, 1.1)

# ── 9.5  Adiabatyczna temperatura spalania T_ad(λ) — Cantera ──
print("\n  Liczenie T_ad i S_L dla zakresu λ (Cantera)... to może chwilę potrwać.")
lam_arr = np.arange(0.7, 1.81, 0.1)
Tad_arr = np.full_like(lam_arr, np.nan)
SL_arr  = np.full_like(lam_arr, np.nan)
for i, lam_i in enumerate(lam_arr):
    try:
        r = analiza_spalania_cantera(lam=lam_i, T1=T1, p1=p1)
        Tad_arr[i] = r['T_ad']
        SL_arr[i]  = r['S_L'] * 100   # cm/s
    except ct.CanteraError:
        pass

ax5 = fig.add_subplot(gs[2, 0])
ax5.plot(lam_arr, Tad_arr, 'o-', color='#b8463f', linewidth=2, markersize=4)
ax5.axvline(lam, color='gray', linestyle='--', linewidth=1, label=f'λ = {lam} (przyjęte)')
ax5.set_xlabel('Współczynnik nadmiaru powietrza λ [-]', fontsize=11)
ax5.set_ylabel('T_ad [K]', fontsize=11)
ax5.set_title(f'Adiabatyczna temp. spalania T_ad(λ) — {PALIWO_NAZWA}/powietrze', fontsize=10)
ax5.legend(fontsize=9)
ax5.grid(True, alpha=0.3)

# ── 9.6  Prędkość spalania laminarnego S_L(λ) — Cantera ──
ax6 = fig.add_subplot(gs[2, 1])
ax6.plot(lam_arr, SL_arr, 'o-', color='#2f7a3f', linewidth=2, markersize=4)
ax6.axvline(lam, color='gray', linestyle='--', linewidth=1, label=f'λ = {lam} (przyjęte)')
ax6.set_xlabel('Współczynnik nadmiaru powietrza λ [-]', fontsize=11)
ax6.set_ylabel('S_L [cm/s]', fontsize=11)
ax6.set_title(f'Prędkość spalania laminarnego S_L(λ) — {PALIWO_NAZWA}/powietrze', fontsize=10)
ax6.legend(fontsize=9)
ax6.grid(True, alpha=0.3)

plt.savefig('wyniki_silnik.png', dpi=150, bbox_inches='tight')
print("\n  Wykres zapisany: wyniki_silnik.png")
plt.show()

# ─────────────────────────────────────────────────────────────
# 10. ANALIZA PARAMETRYCZNA — wpływ ε na IMEP i η_th
# ─────────────────────────────────────────────────────────────

print("\n  Analiza parametryczna — wpływ stopnia sprężania ε:")
print(f"  {'ε':>5}  {'p_max [bar]':>12}  {'IMEP [bar]':>11}  {'η_th [%]':>9}")
print("  " + "-" * 44)

eps_arr = np.arange(6, 23, 1) if TYP_SILNIKA == 'ZS' else np.arange(6, 15, 1)
for e in eps_arr:
    Vc_e = Vd / (e - 1)
    w = symuluj_obieg(
        eps=e, p1=p1, T1=T1, lam=lam,
        m=m, phi0=phi0, delta_phi=delta_phi, eta_i=eta_i,
        Hu=Hu, L_st=L_st, Vd=Vd, Vc=Vc_e,
        gamma=gamma, R_air=R_air, r_kl=r_kl, dphi=1.0
    )
    print(f"  {e:>5.0f}  {w['p_max']/1e5:>12.2f}  {w['IMEP']/1e5:>11.3f}  {w['eta_th']*100:>9.2f}")