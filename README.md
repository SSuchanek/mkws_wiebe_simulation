# Symulacja termodynamiczna silnika tłokowego — model Wiebescha

Projekt na przedmiot **Metody Komputerowe w Spalaniu** (Politechnika Warszawska).

Numeryczna symulacja obiegu pracy silnika tłokowego ZI (Otto) i ZS (Diesel) z modelem spalania Wiebescha. Program całkuje równania termodynamiczne krok po kroku po kącie obrotu wału korbowego i generuje wykresy oraz wskaźniki pracy silnika.

---

## Fizyka modelu

### Geometria układu korbowo-tłokowego

Objętość cylindra w funkcji kąta OWK wyznaczana jest z pełnej kinematyki korbowodu:

$$V(\varphi) = V_c + \frac{V_s}{2} \cdot \frac{1 - \cos\varphi + \frac{1}{\lambda_r}\left(1 - \sqrt{1 - \lambda_r^2 \sin^2\varphi}\right)}{1 + \frac{1}{\lambda_r}}$$

gdzie $\lambda_r = r/l$ to stosunek ramienia korby do długości korbowodu.

### Model spalania Wiebescha

Udział masy spalonego paliwa opisuje funkcja:

$$x(\varphi) = 1 - \exp\!\left(-a \cdot \left(\frac{\varphi - \varphi_0}{\Delta\varphi}\right)^{m+1}\right)$$

- $a = 6.908$ — stała kształtu (odpowiada $\eta_v = 99.9\%$)
- $\varphi_0$ — kąt początku spalania [°OWK]
- $\Delta\varphi$ — czas trwania spalania [°OWK]
- $m$ — wykładnik kształtu (ZI: $m \approx 2$, ZS: $m \approx 0.8$)

### Bilans energii (I zasada termodynamiki)

Ciśnienie całkowane jest krokowo metodą różnic skończonych:

$$\frac{dp}{d\varphi} = \frac{\gamma - 1}{V} \frac{dQ}{d\varphi} - \frac{\gamma \, p}{V} \frac{dV}{d\varphi}$$

Przyrost ciepła chemicznego w kroku $d\varphi$:

$$dQ = \frac{dx}{d\varphi} \cdot Q_{total} \cdot d\varphi, \qquad Q_{total} = m_f \cdot H_u \cdot \eta_i$$

---

## Wyniki i wskaźniki

| Wskaźnik | Opis |
|---|---|
| $p_{max}$ | Maksymalne ciśnienie w cylindrze [bar] |
| $T_{max}$ | Maksymalna temperatura gazów [K] |
| IMEP | Średnie ciśnienie indykowane [bar] |
| $\eta_{th}$ | Sprawność termodynamiczna [%] |

Program generuje cztery wykresy:

- **p-V** — wykres indykatorowy (pole = praca indykowana)
- **p(φ)** — ciśnienie w funkcji kąta OWK
- **T(φ)** — temperatura gazów w funkcji kąta OWK
- **x(φ)** — funkcja Wiebescha dla różnych wykładników $m$

Oraz tabelę analizy parametrycznej — wpływ stopnia sprężania $\varepsilon$ na IMEP i $\eta_{th}$.

---

## Instalacja i uruchomienie

```bash
pip install numpy matplotlib
python wiebe_engine_simulation.py
```

Wymagany Python ≥ 3.10 i NumPy ≥ 2.0.

---

## Parametry do zmiany

Wszystkie parametry znajdują się na początku pliku w sekcjach 1–2:

```python
TYP_SILNIKA = 'ZI'   # 'ZI' lub 'ZS'
eps         = 10.0   # stopień sprężania
p1          = 1.01325e5  # ciśnienie dolotowe [Pa]
T1          = 320.0  # temperatura dolotowa [K]
lam         = 1.0    # współczynnik nadmiaru powietrza
m           = 2.0    # wykładnik kształtu Wiebescha
phi0        = -10.0  # kąt początku spalania [°OWK]
delta_phi   = 60.0   # czas trwania spalania [°OWK]
eta_i       = 0.95   # sprawność wskazana
```

---

## Struktura pliku

```
wiebe_engine_simulation.py
│
├── 1. Parametry silnika
├── 2. Parametry modelu Wiebescha
├── 3. Stałe fizyczne i dane paliwa
├── 4. Geometria układu korbowo-tłokowego
│     └── objetosc_dokladna(phi_deg)
├── 5. Funkcja Wiebescha
│     └── wiebe(phi_deg, phi0, delta_phi, m)
├── 6. Symulacja obiegu — całkowanie RK/FD
│     └── symuluj_obieg(...)
├── 7. Uruchomienie symulacji
├── 8. Wydruk wskaźników
├── 9. Wykresy (matplotlib)
└── 10. Analiza parametryczna (wpływ ε)
```

---

## Literatura

- Rychter T., Teodorczyk A., *Teoria silników tłokowych*, WKŁ, Warszawa 2006
- Heywood J.B., *Internal Combustion Engine Fundamentals*, McGraw-Hill, 1988
- Wiebe I.I., *Brennverlauf und Kreisprozess von Verbrennungsmotoren*, VEB Verlag Technik, Berlin 1970
