# 03. Wykryj wspólne pochodzenie danych w podziale

Status: ready-for-agent

State: open

Blocked by: brak

## Co dostarcza karta

Badacz sprawdza przypisanie przypadków do zbiorów przed porównaniem algorytmów. Otrzymuje listę grup, wariantów i wykrytych konfliktów.

## Kryteria odbioru

- [ ] Warianty jednego katastrowego ciągu, fragmenty tego samego systemu, powtórne skany i warianty analitycznej geometrii mają wspólną grupę pochodzenia.
- [ ] Próba przydzielenia dwóch wariantów jednej grupy do różnych zbiorów wykazuje konflikt z identyfikatorami obu przypadków. Spójny podział odtwarza się z zapisanego opisu.
- [ ] Dwa oglądane LAZ, fragment Kalackiej i znane przypadki syntetyczne są oznaczone jako materiał rozwojowy. Raport jawnie pokazuje, gdy niezależny zbiór testowy jest pusty.
- [ ] Identyfikacja przypadków i wersji podziału działa lokalnie bez kopiowania całych chmur. Brak potwierdzenia, że dwa skany reprezentują różne jaskinie, jest widoczny w opisie.
