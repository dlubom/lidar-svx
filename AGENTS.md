# lidar-svx

Biblioteka wyznacza z chmury LiDAR sieć korytarzy i zapisuje ją w Survex. Przy zmianie algorytmu oceniaj wynik na tych samych przypadkach co poprzednią wersję. Wersja źródła, jednostki i wzorzec muszą być możliwe do odtworzenia.

## Agent skills

### Issue tracker

Specyfikacje i zadania prowadzimy lokalnie w `.scratch/`. Przy planowaniu użyj `to-spec` i `to-tickets`. Przy realizacji odczytaj wskazaną kartę, jej zależności i specyfikację nadrzędną. [Konwencja trackera](docs/agents/issue-tracker.md)

### Triage labels

Przy kwalifikowaniu zadania stosuj pięć domyślnych etykiet Matta. Etykieta określa, kto może podjąć pracę. [Znaczenie etykiet](docs/agents/triage-labels.md)

### Domain docs

Repozytorium ma jeden [słownik dziedziny](GLOSSARY.md). Przy zmianie sposobu rekonstrukcji, metryk lub podziału danych przeczytaj [metodykę oceny szkieletu](docs/badania/ocena-szkieletu.md). Źródła danych, próbę na Kalackiej i zasady korzystania z dokumentacji opisują [instrukcje dziedzinowe](docs/agents/domain.md).

## Weryfikacja

Przy zmianie eksportu skompiluj wynik programem `cavern` i odczytaj stacje, odcinki oraz ewentualne przekroje przez `dump3d`. Przy zmianie ekstrakcji lub metryk pokaż również wynik geometryczny na danych o znanym wzorcu. [Polecenia kontroli jakości](README.md#kontrole-jakości)
