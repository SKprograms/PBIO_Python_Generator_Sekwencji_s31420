# s31420, 2026-05-08
# Ten kod to generator losowych sekwencji nukleotydowych z zapisem w formacie FASTA
# Jest to wykonywane dzięki podaną przez użytkownika długością sekwencji, identyfikatorem oraz opisem
# Obliczane są statystki takie jak procentowy udział danego nukleotydu oraz zawartość GC
# Dodatkowymi rzeczami, które są wykonywane są:
#   Walidacja procentowego udziału nukleotydów
#   Transkrypcja 'in silico' wraz z zapisem wygenerowanej sekwencji mRNA w pliku FASTA
#   Translacja wygenerowanej sekwencji na sekwencje aminokwasową (na podstawie tablicy kodonów)
#   Walidacja pliku FASTA użytkownika i sprawdzenie formatu (nagłowek, znaki, szerokość lini) i raportem błędów

import random
import os


# Wymagania podstawowe --------------------------------------------------

# Generowanie losowej sekwencji DNA o długości użytkownika wraz z dodatkowym podanym rozkładem nukleotydów
def generate_sequence(length: int,
                      wagi: dict = None) -> str:
    """Zwraca losową sekwencję DNA o zadanej długości wraz z możliwie podanym rozkładem nukleotydów"""
    nukleotydy = ['A', 'C', 'G', 'T']

    # Dodatkowa funkcjonalność (rozkład nukleotydów od użytkownika)
    if wagi != None:
        lista_wag = [wagi['A'], wagi['C'], wagi['G'], wagi['T']]
        lista_sekwencji = random.choices(nukleotydy, weights=lista_wag, k=length)
    else:
        lista_sekwencji = []
        for i in range(length):
            lista_sekwencji.append(random.choice(nukleotydy))

    wynik_sekwencja = ""
    for znak in lista_sekwencji:
        wynik_sekwencja = wynik_sekwencja + znak

    return wynik_sekwencja


# Obliczanie statystyk wygenerowanej sekwencji DNA
def calculate_stats(sequence: str) -> dict:
    """Zwraca słownik ze statystykami sekwencji dla poszczególnych nukleotydów oraz GC"""
    sekw = sequence.upper()
    dlugosc = len(sekw)

    if dlugosc == 0:
        return {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0, "GC": 0.0}

    statystyki = {
        "A": (sekw.count("A") / dlugosc) * 100,
        "C": (sekw.count("C") / dlugosc) * 100,
        "G": (sekw.count("G") / dlugosc) * 100,
        "T": (sekw.count("T") / dlugosc) * 100,
    }
    statystyki["GC"] = statystyki["G"] + statystyki["C"]

    return statystyki


# Wstawienie podanego imienia w losowe miejsce wygenerowanej sekwencji
def insert_name(sequence: str, name: str) -> str:
    """Wstawia imię w losową pozycję sekwencji zapisane małymi literami"""
    if name == "":
        return sequence
    pozycja = random.randint(0, len(sequence))
    return sequence[:pozycja] + name.lower() + sequence[pozycja:]


# Sformatowanie ciągu znaków do formatu zapisu pliku FASTA
def format_fasta(seq_id: str, description: str,
                 sequence: str, line_width: int = 80) -> str:
    """Zwraca sformatowany rekord FASTA jako string z łamaniem linii co line_width znaków"""
    naglowek = (">" + seq_id + " " + description).strip()

    linie = [naglowek]
    for i in range(0, len(sequence), line_width):
        linie.append(sequence[i:i + line_width])

    wynik_tekst = ""
    for linia in linie:
        if wynik_tekst == "":
            wynik_tekst = linia
        else:
            wynik_tekst = wynik_tekst + "\n" + linia

    return wynik_tekst


# Pobieranie długości sekwencji oraz walidacja podanej liczby
def validate_positive_int(prompt: str,
                          min_val: int = 1,
                          max_val: int = 100_000) -> int:
    """Pobiera od użytkownika liczbę całkowitą z zakresu i powtarza w przypadku błędu"""
    while True:
        try:
            wartosc = int(input(prompt))
            if min_val <= wartosc <= max_val:
                return wartosc
            else:
                print("Wystąpił błąd ! -> wartość musi być liczbą całkowitą z zakresu [" + str(min_val) + ", " + str(
                    max_val) + "]")
        except ValueError:
            print("Wystąpił błąd ! -> wartość musi być liczbą całkowitą z zakresu [" + str(min_val) + ", " + str(
                max_val) + "]")


# Funkcje dodatkowe --------------------------------------------------

# Pobranie od użytkownika procentowego udziału każdego nukleotydu
def pobierz_rozklad_nukleotydow() -> dict | None:
    """Konfigurowalny rozkład nukleotydów sprawdzający czy ich suma wynosi 100%"""
    print("\n")
    print("--- Konfiguracja rozkładu nukleotydów (Opcjonalnie) ---")
    wybor = input("Czy chcesz podać własny procentowy rozkład nukleotydów? (tak/nie): ").strip().lower()

    if wybor != 'tak':
        return None

    while True:
        try:
            print("Podaj wartości procentowe (np. 25). Suma musi wynosić 100")
            a = float(input(" % A: "))
            c = float(input(" % C: "))
            g = float(input(" % G: "))
            t = float(input(" % T: "))

            if a < 0 or c < 0 or g < 0 or t < 0:
                print("Wystąpił błąd ! -> Wartości procentowe nie mogą być ujemne. Spróbuj ponownie")
                continue

            if abs((a + c + g + t) - 100.0) < 0.001:
                return {'A': a, 'C': c, 'G': g, 'T': t}
            else:
                print("Wystąpił błąd ! -> Suma wynosi " + str(
                    a + c + g + t) + "%, a powinna wynosić 100%. Spróbuj ponownie")
        except ValueError:
            print("Wystąpił błąd ! -> Należy podać wartości liczbowe")


# Transkrypcja sekwencji in silico
def transkrypcja(sekwencja: str) -> str:
    """Transkrypcja in silico zamieniająca tyminę na uracyl"""
    return sekwencja.replace('T', 'U').replace('t', 'u')


# Translacja wygenerowanej sekwencji DNA na sekwencję aminokwasową
def translacja_dna_na_bialko(sekwencja: str) -> str:
    """Translacja przepisująca sekwencję DNA na sekwencję aminokwasową wg tablicy kodonów"""
    tabela_kodonow = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                      'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                      'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                      'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                      'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                      'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                      'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                      'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}

    bialko = []
    sekw = sekwencja.upper()
    for i in range(0, len(sekw) - len(sekw) % 3, 3):
        kodon = sekw[i:i + 3]
        bialko.append(tabela_kodonow.get(kodon, '?'))

    wynik_bialko = ""
    for aminokwas in bialko:
        wynik_bialko = wynik_bialko + aminokwas

    return wynik_bialko


# Walidacja istniejącego pliku FASTA pod względem poprawnego formatu
def walidacja_fasta(sciezka_pliku: str):
    """Walidator pliku FASTA sprawdzający nagłówek, długość linii i znaki"""
    print("\n")
    print("--- Rozpoczynam walidację pliku " + sciezka_pliku + " ---")

    if os.path.exists(sciezka_pliku) == False:
        print("Wystąpił błąd ! -> Podany plik nie istnieje")
        return

    with open(sciezka_pliku, 'r') as plik:
        linie = plik.readlines()

    if len(linie) == 0:
        print("Plik jest pusty")
        return

    bledy = []

    if linie[0].startswith('>') == False:
        bledy.append("Brak poprawnego nagłówka (pierwsza linia nie zaczyna się od '>')")

    w_sekwencji = False
    for numer_linii, linia in enumerate(linie):
        linia = linia.strip()

        if linia == "":
            continue

        if linia.startswith('>') == True:
            w_sekwencji = True
            continue

        if w_sekwencji == True:
            if len(linia) > 80:
                bledy.append("Linia " + str(numer_linii + 1) + " przekracza 80 znaków (ma " + str(len(linia)) + ")")

            if linia.replace("*", "").isalpha() == False:
                bledy.append("Linia " + str(numer_linii + 1) + " zawiera niedozwolone znaki (niebędące literami)")

    if len(bledy) > 0:
        print("Znaleziono następujące błędy w pliku:")
        for blad in bledy:
            print(" - " + blad)
    else:
        print("Plik FASTA jest poprawny")


# Funkcja główna --------------------------------------------------

# Wywoływanie funkcjonalności
def main():
    """Główna funkcja programu integrująca wszystkie etapy generowania i przetwarzania sekwencji"""

    print("\n")
    print("--- Generator sekwencji DNA i walidator FASTA ---")
    print("Co chcesz zrobić?")
    print("1. Wygeneruj nową sekwencję DNA (i zapisz do FASTA)")
    print("2. Sprawdź poprawność istniejącego pliku FASTA")

    wybor_start = input("Wybierz opcję (1 lub 2): ").strip()

    if wybor_start == "1":
        # Krok 1: Wstępna konfiguracja
        wagi = pobierz_rozklad_nukleotydow()

        print("\n")
        print("--- Generowanie nowej sekwencji ---")
        dlugosc = validate_positive_int("Podaj długość sekwencji: ", 1, 100_000)

        # Walidacja ID
        while True:
            id_sekwencji = input("Podaj ID sekwencji: ").strip()
            if id_sekwencji == "":
                print("Wystąpił błąd ! -> ID nie może być puste ani składać się z samych spacji")
            elif " " in id_sekwencji or "\t" in id_sekwencji:
                print("Wystąpił błąd ! -> ID sekwencji nie może zawierać białych znaków w środku")
            else:
                break

        # Opis jest opcjonalny
        opis = input("Podaj opis sekwencji (opcjonalne): ").strip()

        # Walidacja imienia
        while True:
            imie = input("Podaj imię: ").strip()
            if imie == "":
                print("Wystąpił błąd ! -> Imię nie może być puste ani składać się z samych spacji")
            elif imie.isalpha() == False:
                print("Wystąpił błąd ! -> Imię może zawierać tylko litery (bez cyfr i znaków specjalnych)")
            else:
                break

        # Krok 2: Generowanie i modyfikacja
        sekwencja_bazowa = generate_sequence(dlugosc, wagi)
        statystyki = calculate_stats(sekwencja_bazowa)
        sekwencja_koncowa = insert_name(sekwencja_bazowa, imie)

        # Krok 3: Formatowanie
        fasta_dna = format_fasta(id_sekwencji, opis, sekwencja_koncowa)
        sekwencja_mrna = transkrypcja(sekwencja_bazowa)
        fasta_mrna = format_fasta(id_sekwencji + "_mRNA", "Transkrypcja in silico", sekwencja_mrna)
        sekwencja_bialka = translacja_dna_na_bialko(sekwencja_bazowa)
        fasta_bialko = format_fasta(id_sekwencji + "_Protein", "Translacja na bialko", sekwencja_bialka)

        # Krok 4: Zapis do pliku
        nazwa_pliku = id_sekwencji + ".fasta"
        with open(nazwa_pliku, "w") as plik:
            plik.write(fasta_dna + "\n\n")
            plik.write(fasta_mrna + "\n\n")
            plik.write(fasta_bialko + "\n")

        print("\n")
        print("Sekwencja zapisana do pliku: " + nazwa_pliku)

        # Krok 5: Wypisanie statystyk
        print("\n")
        print("Statystyki sekwencji (n=" + str(dlugosc) + "):")
        print("  A: " + str(round(statystyki['A'], 2)) + "%")
        print("  C: " + str(round(statystyki['C'], 2)) + "%")
        print("  G: " + str(round(statystyki['G'], 2)) + "%")
        print("  T: " + str(round(statystyki['T'], 2)) + "%")
        print("  GC-content: " + str(round(statystyki['GC'], 2)) + "%")

        # Krok 6: Opcjonalna walidacja po wygenerowaniu
        print("\n")
        print("--- Walidator plików FASTA ---")
        wybor_walidacji = input("Czy chcesz sprawdzić wygenerowany plik walidatorem? (tak/nie): ").strip().lower()
        if wybor_walidacji == 'tak':
            sciezka_do_sprawdzenia = input(
                "Podaj nazwę pliku do sprawdzenia ((wraz z .fasta), ostatni wygenerowany -> " + nazwa_pliku + ")): ").strip()
            if sciezka_do_sprawdzenia == "":
                sciezka_do_sprawdzenia = nazwa_pliku
            walidacja_fasta(sciezka_do_sprawdzenia)

    elif wybor_start == "2":
        # Użytkownik wybrał tylko opcję drugą (walidator)
        print("\n")
        sciezka_do_sprawdzenia = input("Podaj nazwę pliku do sprawdzenia (wraz z .fasta): ").strip()
        if sciezka_do_sprawdzenia == "":
            print("Wystąpił błąd ! -> Nazwa pliku nie może być pusta")
        else:
            walidacja_fasta(sciezka_do_sprawdzenia)

    else:
        print("\n")
        print("Wystąpił błąd ! -> Wybrano niepoprawną opcję. Uruchom program ponownie")


if __name__ == "__main__":
    main()