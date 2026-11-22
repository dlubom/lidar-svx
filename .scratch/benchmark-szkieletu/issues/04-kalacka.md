# 04. Oceń rekonstrukcję Kalackiej po kompilacji

Status: ready-for-agent

State: open

Blocked by: 01, 02, 03

## Co dostarcza karta

Badacz uruchamia pełną próbę na fragmencie Kalackiej i dostaje raport przypisany do wersji danych, generatora, ekstraktora oraz metryk.

## Kryteria odbioru

- [ ] Z rzeczywistych par końców odcinków wyrównanego modelu odtwarzamy trzy wybrane ramiona i syntetyczne XYZ. Przy udokumentowanym ziarnie chmura zgadza się z zapisanym SHA-256.
- [ ] Ekstraktor dostaje tylko chmurę. Po kompilacji Survex sprawdzamy zgodność nazw, położeń i par końców odcinków odczytanych z modelu 3D, z uwzględnieniem ramy lokalnej i precyzji zapisu.
- [ ] Raport jednego przypadku pokazuje grupę danych, wyniki obu metryk, zakres wzorca oraz konkretne braki sieci. Rozróżnia próbkę porównawczą pokrycia osi od dopasowania rzeczywistych przejść.
- [ ] Zapisany wynik można ocenić ponownie tym samym poleceniem oceny bez uruchamiania ekstraktora. Dla zachowanego wyniku odniesienia odczyt daje 18 stacji, 17 odcinków i dwa końce przy trzech końcach wzorca.
