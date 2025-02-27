import math

metadata = {
    'protocolName': 'Cherrypicking and Normalization',
    'author': 'Chaz <protocols@opentrons.com>',
    'source': 'Custom Protocol Request',
    'apiLevel': '2.10'
}


def run(ctx):

    [left_pipette_type, right_pipette_type, d_csv,
     s_csv, diluent_scheme, mix] = get_values(  # noqa: F821
        "left_pipette_type", "right_pipette_type", "d_csv", "s_csv",
        "diluent_scheme", "mix")

    tiprack_map = {
        'p10_single': 'opentrons_96_filtertiprack_10ul',
        'p50_single': 'opentrons_96_filtertiprack_200ul',
        'p300_single_gen1': 'opentrons_96_filtertiprack_200ul',
        'p1000_single_gen1': 'opentrons_96_filtertiprack_1000ul',
        'p20_single_gen2': 'opentrons_96_filtertiprack_20ul',
        'p300_single_gen2': 'opentrons_96_filtertiprack_200ul',
        'p1000_single_gen2': 'opentrons_96_filtertiprack_1000ul'
    }

    # load labware
    transfer_info_d = [[val.strip().lower() for val in line.split(',')]
                       for line in d_csv.splitlines()
                       if line.split(',')[0].strip()][1:]

    transfer_info_s = [[val.strip().lower() for val in line.split(',')]
                       for line in s_csv.splitlines()
                       if line.split(',')[0].strip()][1:]

    for line in transfer_info_d:
        s_lw, s_slot, d_lw, d_slot = line[:2] + line[4:6]
        for slot, lw in zip([s_slot, d_slot], [s_lw, d_lw]):
            if not int(slot) in ctx.loaded_labwares:
                ctx.load_labware(lw.lower(), slot)

    for line in transfer_info_s:
        s_lw, s_slot, d_lw, d_slot = line[:2] + line[4:6]
        for slot, lw in zip([s_slot, d_slot], [s_lw, d_lw]):
            if not int(slot) in ctx.loaded_labwares:
                ctx.load_labware(lw.lower(), slot)

    # load tipracks in remaining slots
    avail_slots = [str(slot) for slot in range(1, 13)
                   if slot not in ctx.loaded_labwares]
    num_avail_slots = len(avail_slots)
    num_pipettes = len([pip for pip in [left_pipette_type, right_pipette_type]
                        if pip])
    if num_pipettes == 0:
        raise Exception('Must select at least 1 pipette.')
    pipettes = {'left': None, 'right': None}
    for i, (pip_type, side) in enumerate(
            zip([left_pipette_type, right_pipette_type], pipettes.keys())):
        if pip_type:
            tiprack_type = tiprack_map[pip_type]
            tipracks = []
            if i == 0:
                num_racks = math.ceil(num_avail_slots/num_pipettes)
                slots = avail_slots[:num_racks]
            else:
                num_racks = math.floor(num_avail_slots/num_pipettes)
                start = num_avail_slots - num_racks
                slots = avail_slots[start:]
            for slot in slots:
                tipracks.append(ctx.load_labware(tiprack_type, str(slot)))
        # load pipette
        pipettes[side] = ctx.load_instrument(pip_type, side,
                                             tip_racks=tipracks)

    tip_log = {
        pip: {'count': 0, 'max': len(pip.tip_racks*96)}
        for pip in pipettes.values() if pip
    }

    def pick_up(pip):
        if tip_log[pip]['count'] == tip_log[pip]['max']:
            ctx.pause('Please refill {}µl tipracks before \
resuming.'.format(pip.max_volume))
            pip.reset_tipracks()
            tip_log[pip]['count'] = 0
        pip.pick_up_tip()
        tip_log[pip]['count'] += 1

    def parse_well(well):
        letter = well[0]
        number = well[1:]
        return letter.upper() + str(int(number))

    ctx.comment('Transferring diluent to wells based on Dilutant CSV...')
    for line in transfer_info_d:
        [_, s_slot, s_well, asp_h, _, d_slot, d_well, disp_h, vol,
         pip] = line[:10]
        pipette = pipettes[pip]
        if diluent_scheme == 'always':
            pick_up(pipette)
        source = ctx.loaded_labwares[
            int(s_slot)].wells_by_name()[parse_well(s_well)]
        dest = ctx.loaded_labwares[
            int(d_slot)].wells_by_name()[parse_well(d_well)]
        if not pipette.has_tip:
            pick_up(pipette)
        pipette.transfer(float(vol), source.bottom(float(asp_h)),
                         dest.bottom(float(disp_h)), new_tip='never')
        if diluent_scheme == 'always':
            pipette.drop_tip()
    for pip in pipettes.values():
        if pip:
            if pip.has_tip:
                pip.drop_tip()

    ctx.comment('Transferring DNA to wells based on Sample CSV...')
    for line in transfer_info_s:
        [_, s_slot, s_well, asp_h, _, d_slot, d_well, disp_h, vol,
         pip] = line[:10]
        source = ctx.loaded_labwares[
            int(s_slot)].wells_by_name()[parse_well(s_well)]
        dest = ctx.loaded_labwares[
            int(d_slot)].wells_by_name()[parse_well(d_well)]
        pipette = pipettes[pip]
        pick_up(pipette)
        pipette.transfer(float(vol), source, dest, new_tip='never')
        if mix:
            pipette.mix(3, float(vol), dest)
        pipette.blow_out()
        pipette.drop_tip()

    ctx.comment('Protocol complete.')
