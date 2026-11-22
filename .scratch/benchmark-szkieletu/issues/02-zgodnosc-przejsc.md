# 02. Policz odtworzone i błędne przejścia

Status: ready-for-agent

State: open

Blocked by: brak

## Co dostarcza karta

Badacz porównuje szkielet z oznaczonymi wylotami i przejściami. Raport wskazuje, gdzie zgubiono odnogę lub dodano przejście.

## Kryteria odbioru

- [ ] Dla wzorca z trzema ramionami i wyniku zachowującego dwa oznaczone przejścia otrzymujemy TP=2, FP=0, FN=1, F1=0,8 oraz identyfikator brakującego ramienia.
- [ ] Gwiazda i dwa rozwidlenia łączące te same cztery wyloty wewnątrz jednej otwartej sali dają zgodny wynik. Dodanie stacji pośredniej również zachowuje wynik.
- [ ] Pominięcie jednego ramienia wokół przegrody i dodanie nadmiarowego obejścia w innym miejscu ujawnia oba błędy, także przy równej liczbie cykli. Jedno przejście wynikowe może pokryć najwyżej jedno wzorcowe.
- [ ] Raport zawiera liczebności, precyzję, czułość, F1 i zakres adnotacji. Brak adnotacji odróżnia od oznaczonego przypadku z pustym wynikiem. Przy zerowym mianowniku wskazuje brak wartości danej miary.
