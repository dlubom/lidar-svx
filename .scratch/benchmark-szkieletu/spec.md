# Benchmark szkieletu jaskini

Status: ready-for-agent

## Problem

Badacz może otrzymać kompilowalny plik Survex, który gubi ramię rozwidlenia albo prowadzi odcinek przez skałę. Odległość od osi pomiarowej zależy od szerokości korytarza i wyboru stacji w sali. Potrzebujemy sposobu porównania kolejnych algorytmów na tych samych jaskiniach i rozpoznania, gdzie wynik wymaga oględzin.

## Wynik dla użytkownika

Dla jednej wersji algorytmu i jednej jaskini benchmark zwraca raport o rzeczywistych połączeniach, przebiegu całych odcinków względem wnętrza i zakresie niezależnego wzorca. Można ponownie ocenić zapisany wynik bez trenowania modelu lub ponownej ekstrakcji. Na tym etapie podajemy wartości bez progu zaliczenia.

## Potrzeby użytkownika

1. Jako badacz chcę znaleźć odcinek wychodzący poza korytarz, nawet jeśli obie jego stacje leżą wewnątrz.
2. Chcę porównać przejścia przez otwartą salę z zachowaniem rzeczywistych wylotów, również wtedy, gdy stacje wypadają w innych miejscach.
3. Chcę wskazać utraconą odnogę, fałszywe połączenie i brakujące obejście wokół filara.
4. Chcę odczytać długość nieocenioną oraz zasięg opisu przejść, gdy skan ma luki lub udokumentowano tylko część jaskini.
5. Chcę porównać metodę geometryczną i model uczony przy jednakowych wejściach, wzorcu i sposobie liczenia miar.
6. Chcę prześledzić wynik od źródłowego pomiaru w `.3d` przez chmurę i `.svx` aż do odcinków ponownie odczytanych przez Survex.
7. Chcę wiedzieć, które niezależne jaskinie zasilały trening, walidację i test oraz czy procedura wykrywa wspólne pochodzenie wariantów.
8. Chcę móc dołączyć przyszły skan i rozpoznać, które miary da się już obliczyć z dostępnych adnotacji.

## Decyzje i sposób odbioru

- Wspólnym punktem odbioru jest uruchomienie benchmarku na zapisanym wyniku i odczyt raportu. Ekstraktor i model ML mają korzystać z tego samego sposobu oceny.
- Wzorzec topologii opisuje realne przejścia po zbudowaniu bryły. Raport liczy TP, FP i FN po dopasowaniu jeden do jednego, podaje precyzję, czułość i F1 oraz wskazuje konkretne pominięte i nadmiarowe przejścia.
- Na całej eksportowanej łamanej mierzymy długość wewnątrz, poza i nieocenioną, z jednostką i zakresem kontroli. Zapisujemy najdłuższe przekroczenie i jego położenie. Maksymalną głębokość można podać, gdy pozwala na to wzorzec.
- Podział obejmuje całe jaskinie, wspólne systemy lub rodziny analitycznych geometrii. Każdy wariant dziedziczy przynależność grupy. Oba już oglądane LAZ i dotychczasowe przypadki syntetyczne są materiałem rozwojowym.
- Zapisujemy identyfikatory danych, wersje algorytmu i sposobu oceny, parametry oraz ziarna symulacji. Brak wzorca, pusty mianownik i brak zwróconej sieci muszą być rozpoznawalne w raporcie.
- Dwie trasy analityczne mają znane odpowiedzi. Trasa w sali daje 0 m poza wnętrzem. Przekątna zakrętu daje `4√2 m` poza korytarzem i sięga 2 m w skałę. Dla topologii kontrolami są rozwidlenie, równoważne połączenia w sali i osobne drogi wokół przegrody.
- Kalacka dostarcza próbę od danych do kompilacji i odczytu przez `dump3d`. Dwa rzeczywiste LAZ sprawdzają zachowanie raportu przy ograniczonych adnotacjach, niepotwierdzonej skali i braku pełnej bryły.

## Karty

1. [Zmierz skrót przez skałę](issues/01-skrot-przez-skale.md)
2. [Policz odtworzone i błędne przejścia](issues/02-zgodnosc-przejsc.md)
3. [Wykryj wspólne pochodzenie danych w podziale](issues/03-podzial-jaskin.md)
4. [Oceń rekonstrukcję Kalackiej po kompilacji](issues/04-kalacka.md)
5. [Pokaż zakres oceny dwóch rzeczywistych LAZ](issues/05-rzeczywiste-laz.md)

## Granica pierwszego etapu

Dobór lub trening nowego ekstraktora, rozpoznawanie początku ciągu i eksport LRUD pozostają tematami kolejnych specyfikacji. Punktem odniesienia pierwszego etapu jest dzisiejszy ekstraktor. [Metodyka i definicje miar](../../docs/badania/ocena-szkieletu.md) · [Źródła i wynik początkowy](../../docs/badania/rozpoznanie.md)
