# Symulacja termodynamiczna silnika tłokowego — model Wiebescha & Cantera

Projekt na przedmiot **Metody Komputerowe w Spalaniu** (Politechnika Warszawska).

Numeryczna symulacja obiegu pracy silnika tłokowego ZI (Otto) i ZS (Diesel) połączona z **analizą kinetyki chemicznej** za pomocą pakietu Cantera. Program całkuje równania termodynamiczne krok po kroku po kącie obrotu wału korbowego (OWK), wyznacza parametry równowagowe oraz strukturę płomienia laminarnego, generując komplet wskaźników i wykresów diagnostycznych.

---

## Fizyka i chemia modelu

### 1. Kinetyka chemiczna i równowaga (Cantera)
Do oceny fizycznej poprawności założeń wykorzystano model paliwa referencyjnego (**CH₄ / metan**) oraz mechanizm reakcji **GRI-3.0** (53 gatunki, 325 reakcji):
* **Adiabatyczna temperatura spalania ($T_{ad}$)** wyznaczana jest przy założeniu stałej enthalpii i ciśnienia (równowaga `HP`) dla zadanego współczynnika nadmiaru powietrza $\lambda$.
* **Równowagowy skład spalin** definiuje profile molowe kluczowych produktów reakcji.
* **Laminarna prędkość spalania ($S_L$)** obliczana jest poprzez rozwiązanie jednowymiarowego płomienia swobodnego (model `FreeFlame`).

### 2. Geometria układu korbowo-tłokowego
Objętość cylindra w funkcji kąta OWK wyznaczana jest z pełnej kinematyki korbowodu:

$$V(\varphi) = V_c + V_d \cdot \frac{(1 - \cos\varphi) + \frac{1}{\lambda_r}\left(1 - \sqrt{1 - \lambda_r^2 \sin^2\varphi}\right)}{2}$$

gdzie $\lambda_r = r/l$ to stosunek ramienia korby do długości korbowodu, $V_c$ to objętość komory spalania, a $V_d$ to pojemność skokowa.

### 3. Model spalania Wiebescha
Udział masy spalonego paliwa opisuje funkcja:

$$x(\varphi) = 1 - \exp\!\left(-a \cdot \left(\frac{\varphi - \varphi_0}{\Delta\varphi}\right)^{m+1}\right)$$

* $a = 6.908$ — stała kształtu (odpowiada wypaleniu mieszanki w $\eta_v = 99.9\%$)
* $\varphi_0$ — kąt początku spalania [°OWK]
* $\Delta\varphi$ — czas trwania (kąt) spalania [°OWK]
* $m$ — wykładnik kształtu (dla ZI: $m \approx 2$, dla ZS: $m \approx 0.8$)

### 4. Bilans energii (I zasada termodynamiki)
Ciśnienie w cylindrze całkowane jest krokowo metodą różnic skończonych dla układu zamkniętego (gaz doskonały):

$$\frac{dp}{d\varphi} = \frac{\gamma - 1}{V} \frac{dQ}{d\varphi} - \frac{\gamma \, p}{V} \frac{dV}{d\varphi}$$

---

## Wyniki i wskaźniki (Konfiguracja bazowa)

Poniższe wyniki zostały wygenerowane dla parametrów: $\varepsilon = 10.0$, $\lambda = 1.0$, $\varphi_0 = -10.0^\circ$ OWK, $\Delta\varphi = 60.0^\circ$ OWK.

| Wskaźnik | Wartość | Jednostka | Opis |
|---|---|---|---|
| $Q$ | 1692.69 | J | Ciepło doprowadzone z paliwa |
| $W$ | 889.00 | J | Praca indykowana (pole powierzchni p-V) |
| $p_{max}$ | 60.59 | bar | Maksymalne ciśnienie w cylindrze |
| $T_{max}$ | 3252.2 | K | Maksymalna temperatura gazów w cylindrze |
| $T_{ad}$ | 2235.3 | K | Adiabatyczna temperatura spalania (Cantera) |
| $S_L$ | 42.18 | cm/s | Laminarna prędkość spalania (Cantera) |
| IMEP | 17.780 | bar | Średnie ciśnienie indykowane |
| $\eta_{th}$ | 52.52 | % | Sprawność termodynamiczna obiegu |

### Równowagowy skład spalin ($\lambda = 1.0$)
Udziały molowe głównych składników spalin (powyżej 0.1%) wyznaczone w Cantera:
* **$N_2$** : 0.7083
* **$H_2O$** : 0.1832
* **$CO_2$** : 0.0850
* **$CO$** : 0.0094
* **$O_2$** : 0.0048
* **$H_2$** : 0.0037
* **$OH$** : 0.0030
* **$NO$** : 0.0020

*Uwaga fizyczna:* Maksymalna temperatura w cylindrze ($T_{max} = 3252.2$ K) jest wyższa niż adiabatyczna temperatura spalania z Cantery ($T_{ad} = 2235.3$ K). Wynika to z faktu, że Cantera oblicza $T_{ad}$ dla mieszanki o parametrach początkowych ($T_1, p_1$), podczas gdy w rzeczywistym silniku spalanie rozpoczyna się w czynniku już silnie sprężonym i podgrzanym przez tłok zbliżający się do GMP.

---

## Wykresy diagnostyczne
Po zakończeniu obliczeń program generuje i zapisuje wielopanelowy arkusz wykresów do pliku `wyniki_silnik.png`. Zawiera on:
1. **Wykres indykatorowy p-V** (pole pod krzywą = praca indykowana).
2. **Przebieg ciśnienia $p(\varphi)$** w funkcji kąta OWK.
3. **Przebieg temperatury $T(\varphi)$** czynnika roboczego.
4. **Funkcję Wiebescha $x(\varphi)$** (porównanie dla różnych wykładników kształtu $m$).
5. **Charakterystykę $T_{ad}(\lambda)$** – wpływ składu mieszanki na temperaturę adiabatyczną.
6. **Charakterystykę $S_L(\lambda)$** – wpływ składu mieszanki na laminarną prędkość płomienia.

---

## Instalacja i uruchomienie

### Wymagania
* Python $\ge$ 3.10
* NumPy (obsługujący składnię `np.trapz` lub `np.trapezoid`)
* Matplotlib
* **Cantera** $\ge$ 3.0

### Instrukcja szybkiego startu
```bash
# Tworzenie środowiska i instalacja pakietów w Anaconda Prompt
conda create --name kws_env python=3.10 numpy matplotlib -y
conda activate kws_env
conda install -c cantera cantera -y

# Uruchomienie symulacji
python projektMKWS.py

Wszystkie parametry znajdują się na początku pliku w sekcjach 1–2:

```python
TYP_SILNIKA = 'ZI'       # 'ZI' (iskrowy, CH4) lub 'ZS' (samoczynny, diesel)
eps         = 10.0       # Stopień sprężania [-]
Vd          = 500e-6     # Pojemność skokowa [m³]
p1          = 1.01325e5  # Ciśnienie na początku sprężania [Pa]
T1          = 320.0      # Temperatura na początku sprężania [K]
lam         = 1.0        # Współczynnik nadmiaru powietrza (λ) [-]

m           = 2.0        # Wykładnik kształtu Wiebescha [-]
phi0        = -10.0      # Kąt początku spalania [°OWK] (przed GMP < 0)
delta_phi   = 60.0       # Czas trwania spalania [°OWK]
eta_i       = 0.95       # Sprawność wskazana (straty ciepła) [-]
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
-   Cantera Object-Oriented Software Tool, https://cantera.org/
