# Od chmury LiDAR do ciągu Survex

`lidar-svx` zamienia zarejestrowaną chmurę powierzchni jaskini w **hipotezę przestrzennej sieci pomiarowej**. Przyjmuje LAS, LAZ oraz tekstowe XYZ. Sam wyznacza końce i rozgałęzienia. Dwie drogi, które rozdzielają się i ponownie łączą, mogą utworzyć w Survex zamkniętą pętlę. Osobne fragmenty dostają osobne dowiązania w lokalnym układzie modelu.

## Uruchomienie

Potrzebujesz Pythona 3.11 lub nowszego i `uv`. Uruchom polecenia z katalogu repozytorium. Podaj ścieżkę do własnego skanu i nowy katalog wyników. Parametr `--voxel` określa bok komórki w jednostkach źródłowego XYZ.

```sh
uv run --frozen lidar-svx skan.laz wyniki/jaskinia --voxel 0.2
```

Otwórz wynik w Aven po skompilowaniu pliku programem `cavern`. Oba narzędzia wchodzą w skład Survex.

```sh
cavern -o wyniki/jaskinia/network.3d wyniki/jaskinia/network.svx
aven wyniki/jaskinia/network.3d
```

Dla pełnego skanu `ALL_0345_ALIGN00.laz` użyto poniższych ustawień. Przy tamtej rozdzielczości program odtworzył przejście przez przewężenie widoczne na przekroju przy `Y = 42`.

```sh
uv run --frozen lidar-svx /ścieżka/do/ALL_0345_ALIGN00.laz wyniki/0345 \
  --voxel 0.12 --min-returns 3 --wall-span 11 --band-width 3.3 --min-branch-length 4 --tolerance 0.3
```

Chmura musi być zarejestrowana we wspólnym kartezjańskim XYZ z jednakową jednostką liniową na wszystkich osiach. Jeśli znasz przelicznik na metry, podaj `--metres-per-unit`. Współrzędne geograficzne w stopniach wymagają wcześniejszej reprojekcji. Dla E57 lub PLY można wyeksportować w CloudCompare oryginalne XYZ jako trzy kolumny rozdzielone białymi znakami. W plikach `.xyz` i `.txt` dozwolone są puste linie i komentarze rozpoczynające się od `#`.

## Jak dobrać parametry

| Parametr | Znaczenie |
|---|---|
| `--voxel` | Bok woksela w jednostkach źródła. Dobierz go tak, żeby w poprzek istotnego korytarza mieściło się kilka komórek, a przeciwległe powierzchnie miały dostatecznie gęste odbicia. Zbyt drobna siatka rozrywa rzadko zeskanowane ściany, zbyt gruba może zatrzeć wąskie przejście. |
| `--min-returns` | Minimalna liczba rzeczywistych rekordów w wokselu ściany. Domyślnie 1. Zwiększenie tłumi izolowane odbicia, ale usuwa również rzadko próbkowane powierzchnie. |
| `--wall-span` | Maksymalny odstęp między przeciwległymi seriami odbić w jednostkach XYZ. Domyślnie 50 wokseli. Wybierz skalę odpowiadającą szerokościom interesujących przejść. Do odtworzenia pierwszej wersji służy `--global-parity`. |
| `--band-width` | Szerokość pasma odległości mierzonej wzdłuż wnętrza. Domyślnie pięć wokseli, minimum jeden. Drobniejsze pasma rozróżniają bliższe rozwidlenia. Szersze pasma upraszczają sale i mogą pochłonąć niewielką pętlę. |
| `--min-branch-length` | Dodatkowe przycinanie krótkich końcówek grafu, mierzone po krokach siatki przed uproszczeniem. Domyślnie 0. |
| `--tolerance` | Dopuszczalne odchylenie osi od uproszczonego odcinka w 3D. Domyślnie 0 zachowuje kroki siatki. Przy dodatniej wartości każda nowa cięciwa musi mieścić się w przyjętym wnętrzu. |
| `--metres-per-unit` | Potwierdzona liczba metrów na jednostkę XYZ. Bez tej opcji Survex wyświetla jedną jednostkę jako jeden nominalny metr, a pola długości fizycznej pozostają puste. |
| `--max-cells` | Limit kandydatów na komórki wnętrza. Długie, kręte składowe są przechowywane jako zbiory zajętych indeksów. Domyślnie 48 mln. Obliczenia wymagają pamięci także na graf i transformację odległości. |

Przy ocenie pętli warto powtórzyć ekstrakcję dla sąsiedniej rozdzielczości i szerokości pasma. Rozmiary dobiera się do gęstości pomiaru oraz najmniejszego interesującego przejścia.

## Jak powstaje graf

1. **Od odbić do siatki.** Program czyta wszystkie rekordy partiami, liczy odbicia w każdym wokselu i zapamiętuje pierwszy rzeczywisty rekord w komórce. Początek siatki, skalę zapisu LAS, offsety, ewentualny CRS i SHA-256 wejścia zapisuje w manifeście. Metodę strumieniowego odczytu dostarcza laspy. [Dokumentacja laspy](https://laspy.readthedocs.io/en/latest/basic.html#chunked-reading)
2. **Od siatki do hipotezy wnętrza.** W każdym z 13 kierunków szukamy serii odbić po obu stronach pustej komórki. Para daje głos, jeśli mieści się w zadanej rozpiętości i liczba serii do co najmniej jednego z końców promienia jest nieparzysta. Dawniej dodatkowe odbicia nawet daleko od badanej komórki mogły unieważnić jej głos. Rdzeń wnętrza wymaga 14 z 27 komórek sąsiedztwa z co najmniej 6 głosami, w tym jednym osiowym. Słabiej podpartą grupę dołączamy tylko wtedy, gdy sąsiaduje ścianami z dwoma osobnymi rdzeniami po co najmniej 8 komórek. Każda dołączona komórka ma co najmniej 3 głosy, w tym jeden osiowy, i 10 sąsiadów w oknie 27 z co najmniej 4 głosami. W rdzeniu pomijamy komórki ścian spełniające próg liczności. W słabiej podpartym połączeniu nawet jedno surowe odbicie wyklucza daną komórkę. Dzięki temu pojedyncze echa w cienkiej przegrodzie pomagają rozdzielić równoległe korytarze. `manifest.json` podaje liczby komórek rdzenia oraz dołączonych komórek.
3. **Od wnętrza do rozgałęzień i pętli.** W każdej objętości program wybiera odległą komórkę startową przez pierwsze przejście po grafie. Od niej oblicza odległości po sześciu sąsiadach, dzieli je na pasma i wydziela spójne regiony każdego pasma. Na granicy jednej pary pasm program odrzuca nadmiarowe kontakty w kolejności od najwęższych. Gdy dwa ramiona biegną osobno przez kolejne pasma i później się spotykają, graf może zachować obie drogi. To dyskretny wariant konstrukcji grafu Reeba. Punkt startowy porządkuje obliczenia. Końce wyjściowe to stacje z jedną krawędzią. [Geodezyjne grafy Reeba w analizie wokseli](https://www.sciencedirect.com/science/article/pii/S0022000017301344)
4. **Od regionów do stacji.** Na styku dwóch regionów program wybiera przejście o największym odstępie od granicy modelowanego wnętrza. W regionie prowadzi do niego drogę po sąsiednich komórkach, premiując krótki przebieg i odstęp od granicy. SciPy wykonuje transformację odległości, wyszukiwanie składowych i algorytm Dijkstry. Głosowanie na wnętrze, budowa pasm, wybór przejść i reguły redukcji są implementacją tego projektu. Krótkie odnogi sięgające od rozwidlenia przez najwyżej dwa pasma podlegają redukcji. [Algorytm Dijkstry w SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csgraph.dijkstra.html)
5. **Od stacji do Survex.** Opcjonalne uproszczenie łańcuchów metodą Ramera–Douglasa–Peuckera w XYZ zachowuje węzły i niezależne cykle. Program zapisuje odcinki jako `*data cartesian from to dx dy dz`. Każda składowa dostaje jedną stację `*fix` we wspólnej ramie modelu. Przesunięcie początku i przelicznik jednostek pozwalają wrócić ze stacji Survex do źródłowego XYZ. [Format danych Survex](https://survex.com/docs/manual/datafile.htm)

## Kontrola na chmurze

W CloudCompare użyj **File → Open**, wczytaj źródłowy LAZ oraz `network.poly` z katalogu wyniku. Oba pliki mają oryginalne współrzędne XYZ. Jeżeli CloudCompare zaproponuje globalne przesunięcie dla dużych współrzędnych, zastosuj to samo dla obu plików. Każdy odcinek ma własną polilinię, więc luki są widoczne.

- `network.svx` zawiera wszystkie składowe. Przy wielu składowych powstaje także `principal.svx` z całą spójną siecią o największej sumie długości odcinków, ze wspólnymi nazwami stacji. Manifest podaje jej numer, udział w długości i liczbę cykli.
- `stations.csv` i `legs.csv` zawierają stacje, oryginalne XYZ, numery składowych oraz rzeczywiste pary połączonych stacji.
- `cycles.csv` podaje dla każdego cyklu bazy numer składowej, zamkniętą sekwencję stacji i obwód. Liczba niezależnych cykli wynosi **odcinki − stacje + składowe**. W jednej sieci można wybrać kilka różnych baz.
- `witnesses.csv` podaje dla każdej stacji pierwszy rekord w wokselu odbić o najbliższym środku. Indeks liczy się od zera, dla XYZ według kolejności wierszy z danymi. Jest to punkt do odszukania w źródle. Dokładny odstęp od wszystkich odbić na całym odcinku sprawdza osobny skrypt.
- `manifest.json` zawiera polecenie, parametry, hasze wejścia i kodu, wersje bibliotek oraz statystyki ekstrakcji. Wygenerowany katalog jest publikowany dopiero po zapisaniu kompletu plików.

Dla LAS i LAZ poniższe polecenie rysuje trzy rzuty i sprawdza odległość od oryginalnych odbić w promieniu 0,25 jednostki od każdego z maksymalnie 100 000 wybranych odcinków. Wynik obejmuje indeks i XYZ najbliższego odbicia znalezionego w promieniu. Dobierz próg do skali danych.

```sh
uv run --frozen python scripts/preview_lidar_svx.py skan.laz \
  wyniki/jaskinia --radius 0.25 --check-legs 100000
```

**Granica interpretacji.** Dwie rozległe, równoległe powierzchnie mogą dać głosy w kilku kierunkach nawet przy braku bocznych ścian. Gęste odbicia od rumoszu mogą tworzyć wyspy i pozorne obejścia. Wąskie przewężenia, luki w obserwacji i dwa połączenia mieszczące się w tej samej parze pasm mogą zniknąć. To, które odnogi i obejścia stanowią osobne korytarze jaskini, wymaga obejrzenia przekrojów albo niezależnej dokumentacji terenowej.

## Samodzielna ocena ciągłości

`scripts/benchmark_continuity.py synthetic` generuje chmury XYZ o znanej geometrii, uruchamia ten sam generator, a następnie porównuje odcinki ze wzorcem. W raporcie sprawdź liczbę połączonych fragmentów i cykli, pokrycie znanej osi przez jeden ciąg oraz odsetek długości odcinków wewnątrz korytarza. Opcja `--holdout` zmienia orientację skanu, próbkowanie i położenie siatki. Argumenty po `--` trafiają do generatora.

```sh
uv run --frozen python scripts/benchmark_continuity.py synthetic --holdout -- \
  --voxel 0.12 --min-returns 3 --wall-span 11 --band-width 3.3 --min-branch-length 4 --tolerance 0.3
```

Dla własnego skanu z dostępną wcześniejszą osią tryb `real` przyjmuje CSV z kolumnami `station,x,y,z,edge_exported`. Dla dwóch wskazanych stacji znajduje ich najbliższe punkty na skończonych odcinkach i sprawdza, czy między nimi prowadzi droga w wyeksportowanej sieci. Błędy położenia końców podaje osobno.

Testy wytwarzają własne dane, także w formacie LAZ. Uruchom je poleceniem `uv run --frozen pytest -q tests`.

Dodatkowa kontrola dwóch wizyt skanera w tym samym przekroju jest dostępna dla plików LAS/LAZ z uporządkowanymi blokami akwizycji i czasem GPS. Opcje pokazują `uv run --frozen python -m lidar_svx.acquisition --help` i `uv run --frozen python scripts/compare_acquisition_visits.py --help`.
