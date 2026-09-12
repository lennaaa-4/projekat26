 Prvi dan: 6.9.2026.

 Paketi potrebni za rad su ucitani u program. Set podataka takodje je ucitan u program i zatim je matrica sa transkriptima definisana kao SeuratObject (anotacija kao metapodaci). Uradjena je kontrola kvaliteta prvo na osnovu procenta mitohondrijalne DNK i onda pomocu napravljenog ViolinPlota u kome su gledani podaci o broju transkripata, broju gena. Filtriranje nije uradjeno jer je iz VlnPlota zakljuceno da su podaci spremni za dalju obradu.

 Drugi dan 7.9.2026.

 Odradjena je logaritamska normalizacija podataka kako bi se stabilizovala varijansa u sirovoj ekspresiji. Odredjeni su najvarijabilniji geni (napravljena je i lista 10 najvarijabilnijih) i na osnovu toga uradjen je VariableFeaturesPlot na kome su oblezeni prvih 10 najvarijabilnijih gena i moja ciljena RNK. VariableFeatures kasnije su skalirani kako bi bile smanjene razliike u ekspresiji gena (tj. dovedene na uporediv nivo). Takodje je odredjena i PCA, analiza kljucnih komponenti, kojom se redukuje dimenzionalnost i odredjene su kljucne bioloske razlike (varijacije) prema kojima su geni podeljeni u PC-eve. Na osnovu PCA napravljen je DimHeatMap za 500 celija u PC_1 komponenti koji sluzi da bi se proverilo kakav je zapravo bioloski signal i da li ima ostataka nekevrste tehnickog suma. Napravljen je ElbowPlot pomocu kog je odredjeno koji ce se PC-evi koristiti u daljoj analizi. Takodje je provereno prisustvo LINC00958 u PC-evima. Nakon toga uradjeno je matematicko klasterovanje celija i uradjeno je pet DimPlotova sa razlicitim rezolucijama koji sluze za podesavanje broja klastera (grupa) koji ce biti stvoreni.

 Treci dan 8.9.2026.

 Nakon analize DimPlotova odluceno je da se pokusa da se podaci vizualizuju na druge nacine jer su prethodno napravljeni plotovi bili prilicno nepregledni. To je uradjeno na tri razlicita nacina: dva puta 3D prikazom, PCA bubble plotom sto podrazumeva vizualizaciju PC_1, PC_2 (na osama) i PC_3 (kroz velicinu tacaka) i 2D plotom na koji je nadodata i treca, z-osa za PC_3. Sprovedena je UMAP analiza i napravljen je UMAP DimPlot radi vizualizacije klastera. Napisan je izvestaj za prva tri dana.

 Cetvrti dan 9.9.2026.

 Uradjen je pseudobulk koji je podrazumevao vizualizaciju, pravljenje i definisanje pocetne matrice zatim pravljenje podgrupa koje su takodje morale biti sredjene pre upotrebe, definisanje i uredjivanje metapodatka; i nakon toga uradjena je DESeq2 analiza, medjutim nakon filtracije i zavrsetka analize celokupnu obradu prosla je samo jedna RNK (LINC00486), te je doneta odluka da se narednih dana to analizira.

 Peti dan 10.9.2026.

 Citanje radova i istrazivanje o LINC00486 koja je jedina prosla analizu diferencijalne ekspresije.

 Sesti dan 11.9.2026. 

 Ponavljanje analize diferencijalne ekspresije. Odluceno je da se dodaju i kontrolne grupe (zdrava bela masa) u analizu te uporeda izgleda kao: chronic active vs control, periplaque vs control a ne periplakna vs hronicna aktivna. Nakon dodavanja kontrola, analiza je uspesna sa 11 diferencijalno eksprimiranih gena izmedju periplakne bele mase i kontrola, dok izmedju hronicno aktivnih lezija i kontrola ima 56 diferencijalno eksprimiranih gena. Za ove analize uradjeni su VolcanoPlotovi radi vizualizacije.

 Sedmi dan 12.9.2026.

 Pisanje izvestaja i komentara za figure. Procena dijagnosticke vrednosti kao biomarkera.
