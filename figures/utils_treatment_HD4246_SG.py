anno_leiden_dict = {
    # 'filename': 'HD4246_merged_harmony_res04',
    'resolution': 0.4,
    'labels': {
        '0': 'revSC', #'CLU+',
        '1': 'Fetal', #TACSTD2+ (fetal-like)',
        '2': 'Prolif. ISC', #'LGR5+ prolif.',
        '3':  'HRC', # 'EMP1+ (non-prolif.)',
        '4': 'Ribosomal',
        '5': 'Prolif.',
        '6': 'ISC', # 'LGR5+ (non-prolif.)',
        '7': 'Fetal/HRC', # 'TACSTD2+ EMP1+ (non-prolif.)',
        '8': 'EMT', # 'Undef.',
    }
}

color_code = {
    'treatment_resp': {
        'order': ['down', 'up', 'stable',],
        'colors': ['#5B84B1FF', '#FC766AFF', 'grey']
    }
}