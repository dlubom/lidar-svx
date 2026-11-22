# Ciąg pomiarowy z chmury LiDAR

**Stacja** to punkt ciągu, w którym może zaczynać się lub kończyć odcinek pomiarowy. Stacje wygenerowane z chmury mają współrzędne oszacowane z geometrii skanu.

**Początek ciągu** to wybrana stacja odniesienia do opisu sieci. Może odpowiadać otworowi jaskini, początkowi skanowania lub innemu miejscu wskazanemu przez użytkownika.

**Otwór jaskini** to miejsce, w którym korytarz łączy się z powierzchnią. Jedna jaskinia może mieć kilka otworów.

**Wylot korytarza do sali** to miejsce przejścia między korytarzem a wnętrzem sali. W modelu o znanej geometrii można oznaczyć je przekrojem.

**Szkielet jaskini** to przestrzenna sieć odcinków odwzorowująca przebieg i połączenia korytarzy. W otwartej sali różne przebiegi szkieletu mogą prowadzić do tych samych wylotów.

**Początek skanowania** to pozycja skanera w chwili rozpoczęcia rejestracji. Współrzędne najwcześniejszego odbicia wskazują punkt trafiony wiązką.

**Luka obserwacji** to obszar, w którym pokrycie skanem nie wystarcza do rozstrzygnięcia, czy dwa fragmenty korytarzy są połączone.

**LRUD** to cztery odległości od stacji do lewej ściany, prawej ściany, stropu i spągu. Określenie lewej i prawej strony wymaga podania kierunku przebiegu korytarza.
