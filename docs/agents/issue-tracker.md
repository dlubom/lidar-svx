# Zadania w lokalnym Markdown

Jedno przedsięwzięcie zajmuje jeden katalog `.scratch/<temat>/`. Specyfikacja to `spec.md`, a poszczególne karty to `issues/NN-<temat>.md`. Numeracja zaczyna się od `01` w obrębie specyfikacji. Publikacja przez `to-spec` albo `to-tickets` oznacza zapis w tym katalogu.

- `to-spec` zbiera problem, oczekiwane zachowanie i sposób sprawdzenia wyniku. `to-tickets` dzieli specyfikację na karty, z których każda dostarcza działający i sprawdzalny rezultat.
- Pole `Blocked by` wskazuje numery kart z tej samej specyfikacji. `Status` zawiera [etykietę triage](triage-labels.md), a `State` ma wartość `open` albo `closed`. Agent może podjąć otwartą kartę z etykietą `ready-for-agent`, gdy wszystkie jej zależności są zamknięte.
- Przy zamykaniu karty zaznacz spełnione kryteria i dopisz w sekcji `Wynik` polecenie lub odsyłacz do dowodu odbioru. Zachowaj nazwę pliku, aby odsyłacze dalej działały.

Pierwsze przedsięwzięcie to [benchmark szkieletu jaskini](../../.scratch/benchmark-szkieletu/spec.md).
