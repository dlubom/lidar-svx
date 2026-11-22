# 05. Pokaż zakres oceny dwóch rzeczywistych LAZ

Status: ready-for-agent

State: open

Blocked by: 01, 02, 03

## Co dostarcza karta

Dla każdego z dwóch znanych skanów badacz otrzymuje osobny raport wskazujący dostępny wynik ekstrakcji i zakres niezależnych adnotacji.

## Kryteria odbioru

- [ ] Raport wiąże skan, zapisany lub wygenerowany wynik i ewentualne adnotacje za pomocą identyfikatora źródła. Raporty obu plików zachowują ich osobne układy i jednostki.
- [ ] Jeśli dla danego skanu brak opisu wylotów lub granicy wnętrza, odpowiednia miara jest nieoceniona. Brakujący wynik ekstraktora jest odróżniony od pustej zwróconej sieci. W raporcie można odczytać zakres i przyczynę każdego braku.
- [ ] Na kontrolnym fragmencie z częściową adnotacją miary obejmują wyłącznie ten zakres. Długość nieoceniona pozostaje widoczna także po dołączeniu rozpoznanego fragmentu.
- [ ] Oba skany są raportowane jako materiał rozwojowy. Procedura obsługuje ich lokalne ścieżki i pozwala uzyskać raport bez dołączania dużych plików LAZ do repozytorium.
