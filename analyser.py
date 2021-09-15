from loader import parse_json, load_resource_dir

CODON_TABLE_DIR = load_resource_dir()['codon_table.json']
CODON_TABLE = parse_json(CODON_TABLE_DIR)

SNCT_CODON_TABLE_DIR = load_resource_dir()['snct_codon_table.json']
SNCT_CODON_TABLE = parse_json(SNCT_CODON_TABLE_DIR)

FIVEEND_TO_THIRDEND = '53'
THIRDEND_TO_FIVEEND = '35'
DNA = True
RNA = False


def search_start_codon_index(nucleic_sequence: str, directionality: str) -> int:
    """
    A function that finds start codon index
    returns -1 when there is no start codon
    returns only the first start index found

    :param nucleic_sequence:
    :param directionality:
    :return start_codon_index:

    ex3) directionality = FIVEEEND_TO_THIRDEND
        start codon = AUG
        012 345 678
        ACG AUG CCA GGC UAC GUA
            ↑
        start index is here (3)


    ex4) directionality = THRIDEND_TO_FIVEEND
        start codon = GUA
        012 345 678 9
        AUG UUG GUA GCU CAU UUU
                    ↑
        start index is here (9)

    """
    start_codon_index = -1

    if directionality == FIVEEND_TO_THIRDEND:
        start_codon = "AUG"  # start codon

        for index in range(0, len(nucleic_sequence) - 2):  # search codon string
            single_codon = nucleic_sequence[index: index + 3]
            if single_codon == start_codon:
                start_codon_index = index
                break

    elif directionality == THIRDEND_TO_FIVEEND:
        start_codon = "GUA"  # reversed start codon

        for index in range(len(nucleic_sequence), 2, -1):  # search codon string
            single_codon = nucleic_sequence[index - 3: index]
            if single_codon == start_codon:
                start_codon_index = index
                break

    return start_codon_index


def search_stop_codon_index(nucleic_sequence: str, start_codon_index: int, directionality: str) -> int:
    """
    A Function that finds stop codon index
    returns -1 when there is no stop codon
    returns only the first stop index found

    :param nucleic_sequence:
    :param start_codon_index:
    :param directionality:
    :return:
    """

    stop_codon_index = -1

    if directionality == FIVEEND_TO_THIRDEND:
        stop_codon1, stop_codon2, stop_codon3 = "UAA", "UAG", "UGA"  # stop codon

        for index in range(start_codon_index + 3, len(nucleic_sequence),
                           3):  # search codon string ( 3 digit move per 1 search )
            single_codon = nucleic_sequence[index: index + 3]

            if single_codon == stop_codon1 or single_codon == stop_codon2 or single_codon == stop_codon3:
                stop_codon_index = index
                break

    elif directionality == THIRDEND_TO_FIVEEND:
        stop_codon1, stop_codon2, stop_codon3 = "AAU", "GAU", "AGU"  # stop codon

        for index in range(start_codon_index, 0, -3):  # search codon string
            single_codon = nucleic_sequence[index - 3: index]
            if single_codon == stop_codon1 or single_codon == stop_codon2 or single_codon == stop_codon3:
                stop_codon_index = index
                break

    return stop_codon_index


def is_readable_codon(nucleic_sequence: str) -> bool:
    """
    A function that determines if it is a readable codon

    :param nucleic_sequence:
    :return is_readable_codon:
    """
    # if length of codon string is multiple of 3, it's readable
    if len(nucleic_sequence) % 3 == 0:
        is_readable = True
    else:
        is_readable = False

    return is_readable


def convert_codon_to_amino_acid(codon: str) -> str:
    """
    A function that converts codon into amino acid
    returns '' when there are no such codon
    note: CODON_TABLE is global variable

    :param codon:
    :return amino_acid:
    """
    try:
        amino_acid = CODON_TABLE[codon]
    except KeyError:
        amino_acid = ''

    return amino_acid


def translate_nucleic_sequence(nucleic_sequence: str) -> str:
    """
    A function that converts codon string into amino acid sequence string
    returns "" when all codon in codon string can't be translated

    :param nucleic_sequence:
    :return amino_acid_sequence:
    """
    amino_acid_sequence = ""
    for index in range(0, len(nucleic_sequence), 3):
        codon = nucleic_sequence[index: index + 3]
        if codon:
            amino_acid_sequence += convert_codon_to_amino_acid(codon)

    return amino_acid_sequence


def convert_complementary_base_char(base_char: str, is_dna: bool) -> str:
    """
    A function that converts base char into complementary base char
    returns '' when inputed base doesn't have conplementary base
    :param base_char:
    :param is_dna:
    :return cpm_char:
    """
    dict_of_dna_cpm_base = {'A': 'T',
                            'T': 'A',
                            'C': 'G',
                            'G': 'C'}

    dict_of_rna_cpm_base = {'A': 'U',
                            'U': 'A',
                            'C': 'G',
                            'G': 'C'}

    try:
        if is_dna:
            cpm_char = dict_of_dna_cpm_base[base_char]
        else:
            cpm_char = dict_of_rna_cpm_base[base_char]
    except KeyError:
        cpm_char = ''

    return cpm_char


def convert_complement_base_string(nucleic_sequence: str, is_dna: bool) -> str:
    """
    A function that converts base string into complement base string
    returns "" when all base in base string doesn't have complementary base
    :param nucleic_sequence:
    :param is_dna:
    :return cpm_string:
    """
    cpm_string = ""
    for char in nucleic_sequence:
        cpm_string += convert_complementary_base_char(char, is_dna)

    return cpm_string


def get_complement_nucleic_acid_strand(nucleic_sequence: str, is_dna: bool) -> str:
    """
    A function that returns complement nucleic acid strand
    returns "" when all base in base string doesn't have complementary base
    :param nucleic_sequence:
    :param is_dna:
    :return:
    """
    complement_base_string = convert_complement_base_string(nucleic_sequence, is_dna)
    complement_nucleic_acid_strand = complement_base_string[::-1]

    return complement_nucleic_acid_strand


def transcript_dna_char(base_char: str) -> str:
    """
    A function that transcripts DNA char into mRNA char
    returns '' when base char is not one of DNA bases
    :param base_char:
    :return rna_char:
    """
    dict_of_dna_to_rna_base = {'A': 'U',
                               'T': 'A',
                               'C': 'G',
                               'G': 'C'}
    try:
        rna_char = dict_of_dna_to_rna_base[base_char]

    except KeyError:
        rna_char = ''

    return rna_char


def transcript_dna_string(nucleic_sequence: str) -> str:
    """
    A function that transcripts DNA string info mRNA string
    returns "" when all bases in base string is not one of DNA bases
    :param nucleic_sequence:
    :return rna_string:
    """
    rna_string = ""
    for char in nucleic_sequence:
        rna_string += transcript_dna_char(char)

    return rna_string


def reverse_transcript_rna_char(base_char: str) -> str:
    """
    A function that reverse transcripts mRNA char into DNA char
    returns '' when base char is not one of RNA bases
    :param base_char:
    :return rna_char:
    """
    dict_of_rna_to_dna_base = {'A': 'T',
                               'U': 'A',
                               'C': 'G',
                               'G': 'C'}
    try:
        rna_char = dict_of_rna_to_dna_base[base_char]

    except KeyError:
        rna_char = ''

    return rna_char


def reverse_transcript_rna_string(nucleic_sequence: str) -> str:
    """
    A function that reverse transcripts mRNA string info DNA string
    returns "" when all bases in base string is not one of RNA bases
    :param nucleic_sequence:
    :return rna_string:
    """
    dna_string = ""
    for char in nucleic_sequence:
        dna_string += reverse_transcript_rna_char(char)

    return dna_string


def create_list_of_codon(list_of_base: list) -> list:
    """
    A function that create list of codon with given base chars
    :param list_of_base:
    :return list_of_codon:
    """
    list_of_codon = []
    period = len(list_of_base)
    for i in range(0, period ** 3):  # create all possible codons
        codon = f"{list_of_base[i // period % period]}{list_of_base[i // period % period]}{list_of_base[i % period]}"
        list_of_codon.append(codon)

    return list_of_codon


def search_possible_amino_acids(list_of_base: list) -> list:
    """
    A function that search all possible amino acids with giver base chars
    :param list_of_base:
    :return list of amino_acid:
    """

    set_of_amino_acid = set()
    list_of_codon = create_list_of_codon(list_of_base)
    for codon in list_of_codon:  # convert codon to amino acid
        amino_acid = convert_codon_to_amino_acid(codon)
        set_of_amino_acid.add(amino_acid)  # adding amino acid while deduplicating

    list_of_amino_acid = list(set_of_amino_acid)

    return list_of_amino_acid


def analyze_nucleic_sequence(nucleic_sequence: str, is_dna: bool, directionality: str) -> dict:
    """

    :param nucleic_sequence:
    :param is_dna:
    :param directionality:
    :return dict_of_analyzed_data:
    """

    raw_dna_strand1 = ""
    dna_strand1 = ""
    reversed_dna_strand1 = ""
    raw_dna_strand2 = ""
    dna_strand2 = ""
    reversed_dna_strand2 = ""
    raw_rna_strand1 = ""
    rna_strand1 = ""
    reversed_rna_strand1 = ""
    target_area1 = ""
    raw_translated_string1 = ""
    translated_string1 = ""
    raw_rna_strand2 = ""
    rna_strand2 = ""
    reversed_rna_strand2 = ""
    target_area2 = ""
    raw_translated_string2 = ""
    translated_string2 = ""

    if is_dna:
        if directionality == FIVEEND_TO_THIRDEND:
            raw_dna_strand1 = nucleic_sequence
            dna_strand1 = f"5′ - {raw_dna_strand1} - 3′"
            reversed_dna_strand1 = f"3′ - {raw_dna_strand1[::-1]} - 5′"

            raw_dna_strand2 = get_complement_nucleic_acid_strand(raw_dna_strand1, DNA)
            dna_strand2 = f"3′ - {raw_dna_strand2} - 5′"
            reversed_dna_strand2 = f"5′ - {raw_dna_strand2[::-1]} - 3′"

            raw_rna_strand1 = transcript_dna_string(raw_dna_strand1)
            rna_strand1 = f"3′ - {raw_rna_strand1} - 5′"
            reversed_rna_strand1 = f"5′ - {raw_rna_strand1[::-1]} - 3′"

            raw_rna_strand2 = transcript_dna_string(raw_dna_strand2)
            rna_strand2 = f"5′ - {raw_rna_strand2} - 3′"
            reversed_rna_strand2 = f"3′ - {raw_rna_strand2[::-1]} - 5′"

            start_index1 = search_start_codon_index(raw_rna_strand1, THIRDEND_TO_FIVEEND)
            end_index1 = search_stop_codon_index(raw_rna_strand1, start_index1, THIRDEND_TO_FIVEEND) + 3  # TODO: 계산 다시
            target_area1 = raw_rna_strand1[start_index1: end_index1]
            raw_translated_string1 = translate_nucleic_sequence(target_area1)
            translated_string1 = f" {'  '.join(raw_translated_string1)} "

            start_index2 = search_start_codon_index(raw_rna_strand2, FIVEEND_TO_THIRDEND)
            end_index2 = search_stop_codon_index(raw_rna_strand2, RNA, FIVEEND_TO_THIRDEND)
            target_area2 = raw_rna_strand2[start_index2: end_index2]
            raw_translated_string2 = translate_nucleic_sequence(target_area2)
            translated_string2 = f" {'  '.join(raw_translated_string2)} "

        elif directionality == THIRDEND_TO_FIVEEND:
            raw_dna_strand1 = nucleic_sequence
            dna_strand1 = f"3′ - {raw_dna_strand1} - 5′"  # TODO: 코돈 마다 띄어쓰기 만들기
            reversed_dna_strand1 = f"5′ - {raw_dna_strand1[::-1]} - 3′"

            raw_dna_strand2 = get_complement_nucleic_acid_strand(raw_dna_strand1, DNA)
            dna_strand2 = f"5′ - {raw_dna_strand2} - 3′"
            reversed_dna_strand2 = f"3′ - {raw_dna_strand2[::-1]} - 5′"

            raw_rna_strand1 = transcript_dna_string(raw_dna_strand1)
            rna_strand1 = f"5′ - {raw_rna_strand1} - 3′"
            reversed_rna_strand1 = f"3′ - {raw_rna_strand1[::-1]} - 5′"

            raw_rna_strand2 = transcript_dna_string(raw_dna_strand2)
            rna_strand2 = f"3′ - {raw_rna_strand2} - 3′"
            reversed_rna_strand2 = f"5′ - {raw_rna_strand2[::-1]} - 3′"

            start_index1 = search_start_codon_index(raw_rna_strand1, FIVEEND_TO_THIRDEND)
            end_index1 = search_stop_codon_index(raw_rna_strand1, start_index1, FIVEEND_TO_THIRDEND)
            target_area1 = raw_rna_strand1[start_index1: end_index1]
            raw_translated_string1 = translate_nucleic_sequence(target_area1)
            translated_string1 = f" {'  '.join(raw_translated_string1)} "

            start_index2 = search_start_codon_index(raw_rna_strand2, THIRDEND_TO_FIVEEND)
            end_index2 = search_stop_codon_index(raw_rna_strand2, start_index2, THIRDEND_TO_FIVEEND) + 3
            target_area2 = raw_rna_strand2[start_index2: end_index2]
            raw_translated_string2 = translate_nucleic_sequence(target_area2)
            translated_string2 = f" {'  '.join(raw_translated_string2)} "

    else:
        if directionality == FIVEEND_TO_THIRDEND:
            raw_rna_strand1 = nucleic_sequence
            rna_strand1 = f"5′ - {raw_rna_strand1} - 3′"
            reversed_rna_strand1 = f"3′ - {raw_rna_strand1[::-1]} - 5′"

            raw_dna_strand1 = reverse_transcript_rna_string(raw_rna_strand1)
            dna_strand1 = f"3′ - {raw_dna_strand1} - 5′"
            reversed_dna_strand1 = f"5′ - {raw_dna_strand1[::-1]} - 3′"

            raw_dna_strand2 = convert_complement_base_string(raw_dna_strand1, DNA)
            dna_strand2 = f"5′ - {raw_dna_strand2} - 3′"
            reversed_dna_strand2 = f"3′ - {raw_dna_strand2[::-1]} - 5′"

            raw_rna_strand2 = transcript_dna_string(raw_dna_strand2)
            rna_strand2 = f"3′ - {raw_rna_strand2} - 5′"
            reversed_rna_strand2 = f"5′ - {raw_rna_strand2[::-1]} - 3′"

            start_index1 = search_start_codon_index(raw_rna_strand1, FIVEEND_TO_THIRDEND)
            end_index1 = search_stop_codon_index(raw_rna_strand1, start_index1, FIVEEND_TO_THIRDEND)
            target_area1 = raw_rna_strand1[start_index1: end_index1]
            raw_translated_string1 = translate_nucleic_sequence(target_area1)
            translated_string1 = f" {'  '.join(raw_translated_string1)} "

            start_index2 = search_start_codon_index(raw_rna_strand2, THIRDEND_TO_FIVEEND)
            end_index2 = search_stop_codon_index(raw_rna_strand2, start_index2, THIRDEND_TO_FIVEEND) + 3
            target_area2 = raw_rna_strand2[end_index2 - 3: start_index2]
            raw_translated_string2 = translate_nucleic_sequence(target_area2[::-1])[::-1]
            translated_string2 = f" {'  '.join(raw_translated_string2)} "

        elif directionality == THIRDEND_TO_FIVEEND:
            raw_rna_strand1 = nucleic_sequence
            rna_strand1 = f"3′ - {raw_rna_strand1} - 5′"
            reversed_rna_strand1 = f"5′ - {raw_rna_strand1[::-1]} - 3′"

            raw_dna_strand1 = reverse_transcript_rna_string(raw_rna_strand1)
            dna_strand1 = f"5′ - {raw_dna_strand1} - 3′"
            reversed_dna_strand1 = f"3′ - {raw_dna_strand1[::-1]} - 5′"

            raw_dna_strand2 = convert_complement_base_string(raw_dna_strand1, DNA)
            dna_strand2 = f"3′ - {raw_dna_strand2} - 5′"
            reversed_dna_strand2 = f"5′ - {raw_dna_strand2[::-1]} - 3′"

            raw_rna_strand2 = transcript_dna_string(raw_dna_strand2)
            rna_strand2 = f"5′ - {raw_rna_strand2} - 3′"
            reversed_rna_strand2 = f"3′ - {raw_rna_strand2[::-1]} - 5′"

            start_index1 = search_start_codon_index(raw_rna_strand1, THIRDEND_TO_FIVEEND)
            end_index1 = search_stop_codon_index(raw_rna_strand1, start_index1, THIRDEND_TO_FIVEEND) + 3
            target_area1 = raw_rna_strand1[end_index1 - 3: start_index1]
            raw_translated_string1 = translate_nucleic_sequence(target_area1[::-1])[::-1]
            translated_string1 = f" {'  '.join(raw_translated_string1)} "

            start_index2 = search_start_codon_index(raw_rna_strand2, FIVEEND_TO_THIRDEND)
            end_index2 = search_stop_codon_index(raw_rna_strand2, start_index2, FIVEEND_TO_THIRDEND)
            target_area2 = raw_rna_strand2[start_index2: end_index2]
            raw_translated_string2 = translate_nucleic_sequence(target_area2)
            translated_string2 = f" {'  '.join(raw_translated_string2)} "

    indented_target_area1 = ""
    indented_translated_string1 = ""
    for index in range(0, len(raw_rna_strand1) - len(target_area1)):
        if raw_rna_strand1[index: index + len(target_area1)] == target_area1:
            indented_raw_target_area1 = target_area1
            indented_target_area1 = f"{rna_strand1[:5]}{indented_raw_target_area1}{rna_strand1[-5:]}"
            indented_translated_string1 = f"{rna_strand1[:5]}{translated_string1}{rna_strand1[-5:]}"
            break

    indented_target_area2 = ""
    indented_translated_string2 = ""
    for index in range(0, len(raw_rna_strand2) - len(target_area2)):
        if raw_rna_strand2[index: index + len(target_area2)] == target_area2:
            indented_raw_target_area2 = target_area2
            indented_target_area2 = f"{rna_strand2[:5]}{indented_raw_target_area2}{rna_strand2[-5:]}"
            indented_translated_string2 = f"{rna_strand2[:5]}{translated_string2}{rna_strand2[-5:]}"
            break

    dict_of_analyzed_data = {
        'raw_dna_strand1': raw_dna_strand1,
        'dna_strand1': dna_strand1,
        'reversed_dna_strand1': reversed_dna_strand1,
        'raw_dna_strand2': raw_dna_strand2,
        'dna_strand2': dna_strand2,
        'reversed_dna_strand2': reversed_dna_strand2,
        'raw_rna_strand1': raw_rna_strand1,
        'rna_strand1': rna_strand1,
        'reversed_rna_strand1': reversed_rna_strand1,
        'target_area1': target_area1,
        'indented_target_area1': indented_target_area1,
        'raw_translated_string1': raw_translated_string1,
        'translated_string1': translated_string1,
        'indented_translated_string1': indented_translated_string1,
        'raw_rna_strand2': raw_rna_strand2,
        'rna_strand2': rna_strand2,
        'reversed_rna_strand2': reversed_rna_strand2,
        'target_area2': target_area2,
        'indented_target_area2': indented_target_area2,
        'raw_translated_string2': raw_translated_string2,
        'translated_string2': translated_string2,
        'indented_translated_string2': indented_translated_string2
    }

    return dict_of_analyzed_data


def match_digit_base_to_codon_and_amino_acid(list_of_first_base_boolean: list[bool, bool, bool, bool],
                                             list_of_second_base_boolean: list[bool, bool, bool, bool],
                                             list_of_third_base_boolean: list[bool, bool, bool, bool]) -> list[list, list]:
    """

    :param list_of_first_base_boolean:
    :param list_of_second_base_boolean:
    :param list_of_third_base_boolean:
    :return:
    """
    list_of_selected_first_base = []
    list_of_selected_second_base = []
    list_of_selected_third_base = []
    list_of_base = ['U', 'C', 'A', 'G']
    list_of_matched_codon = []
    list_of_selected_amino_acid = []

    for index in range(4):
        if list_of_first_base_boolean[index]:
            list_of_selected_first_base.append(list_of_base[index])

        if list_of_second_base_boolean[index]:
            list_of_selected_second_base.append(list_of_base[index])

        if list_of_third_base_boolean[index]:
            list_of_selected_third_base.append(list_of_base[index])

    if list_of_selected_first_base and list_of_selected_second_base and list_of_selected_third_base:
        list_of_matched_codon = combine_bases_to_codon(list_of_selected_first_base, list_of_selected_second_base, list_of_selected_third_base)
        list_of_selected_amino_acid = match_codon_with_amino_acid(list_of_matched_codon)

    list_of_codon_and_amino_acid = [list_of_matched_codon, list_of_selected_amino_acid]

    return list_of_codon_and_amino_acid


def combine_bases_to_codon(list_of_selected_first_base: list, list_of_selected_second_base: list,
                           list_of_selected_third_base: list) -> list:
    """

    :param list_of_selected_first_base:
    :param list_of_selected_second_base:
    :param list_of_selected_third_base:
    :return list_of_matched_codon:
    """
    list_of_matched_codon = []

    for first_base in list_of_selected_first_base:
        for second_base in list_of_selected_second_base:
            for third_base in list_of_selected_third_base:
                codon = f"{first_base}{second_base}{third_base}"
                list_of_matched_codon.append(codon)

    return list_of_matched_codon


def match_codon_with_amino_acid(list_of_codon: list) -> list:
    """

    :param list_of_codon:
    :return list_of_matched_amino_acid:
    """
    set_of_matched_amino_acid = set([])
    for codon in list_of_codon:
        matched_amino_acid = SNCT_CODON_TABLE[codon]
        set_of_matched_amino_acid.add(matched_amino_acid)

    list_of_matched_amino_acid = list(set_of_matched_amino_acid)

    return list_of_matched_amino_acid


def combine_all_amino_acids_from_repeated_nucleic_sequence(nucleic_sequence: str, directionality: str) -> list:
    """

    :param nucleic_sequence:
    :param directionality:
    :return list_of_deduplicated_analysis_pair:
    """

    extended_necleic_sequence = f"{nucleic_sequence[-2:]}{nucleic_sequence}{nucleic_sequence[:2]}"
    list_of_codon = []
    list_of_amino_acid = []

    if directionality == FIVEEND_TO_THIRDEND:
        for index in range(0, len(extended_necleic_sequence) - 2):
            codon = extended_necleic_sequence[index: index + 3]
            amino_acid = CODON_TABLE[codon]
            list_of_codon.append(codon)
            list_of_amino_acid.append(amino_acid)

    elif directionality == THIRDEND_TO_FIVEEND:
        for index in range(len(extended_necleic_sequence), 2, -1):
            codon = extended_necleic_sequence[index - 3: index]
            amino_acid = CODON_TABLE[codon[::-1]]
            list_of_codon.append(codon)
            list_of_amino_acid.append(amino_acid)

    list_of_analysis_pair = list(zip(list_of_codon, list_of_amino_acid))

    list_of_deduplicated_analysis_pair = []
    for analysis_pair_ in list_of_analysis_pair:
        if analysis_pair_ not in list_of_deduplicated_analysis_pair:
            list_of_deduplicated_analysis_pair.append(analysis_pair_)

    return list_of_deduplicated_analysis_pair
