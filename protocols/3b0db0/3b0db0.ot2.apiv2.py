import math

metadata = {
    'protocolName': 'Temperature Controlled PCR Prep With Tube Strips',
    'author': 'Rami Farawi <rami.farawi@opentrons.com>',
    'source': 'Custom Protocol Request',
    'apiLevel': '2.7'
}


def run(ctx):

    [num_samp, asp_delay, mix_reps, asp_height,
     p20_mount, m20_mount] = get_values(  # noqa: F821
        "num_samp", "asp_delay", "mix_reps", "asp_height",
        "p20_mount", "m20_mount")

    if not 1 <= num_samp <= 96:
        raise Exception("Enter a sample number between 1-96")

    num_col_floor = math.floor(num_samp/8)
    num_col = math.ceil(num_samp/8)
    remain = num_samp % 8

    # load labware
    temp_mod = ctx.load_module('temperature module gen2', '3')
    sample_plate = temp_mod.load_labware('biorad_96_wellplate_200ul_pcr')
    tip_rack = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot)
                for slot in ['1', '2']]
    saliva = [ctx.load_labware(
            'tuberack_15_tuberack_5000ul', slot)
            for slot in ['5', '6', '7', '8', '9', '10', '11']]
    reagent_plate = ctx.load_labware('pcrstripplate_96_wellplate_200ul', '4')

    # load instrument
    p20 = ctx.load_instrument('p20_single_gen2', p20_mount, tip_racks=tip_rack)
    m20 = ctx.load_instrument('p20_multi_gen2', m20_mount, tip_racks=tip_rack)

    # reagents, setup
    temp_mod.set_temperature(4)
    mastermix = reagent_plate.rows()[0][:2]
    tubes = [tube for tuberack in saliva for tube in tuberack.wells()]

    # add mastermix
    m20.pick_up_tip()
    for i, col in enumerate(sample_plate.rows()[0][:num_col_floor]):
        m20.aspirate(15, mastermix[0] if i < 6 else mastermix[1])
        m20.dispense(15, col)
    m20.drop_tip()

    p20.pick_up_tip()
    for mix, well in zip(reagent_plate.wells()[:remain]
                         if num_col <= 6 else
                         reagent_plate.wells()[8:8+remain],
                         sample_plate.wells()
                         [num_col_floor*8:num_col_floor*8+remain]):
        p20.aspirate(15, mix)
        p20.dispense(15, well)
    p20.drop_tip()

    # add saliva
    airgap = 5
    for s, d in zip(tubes[:num_samp], sample_plate.wells()):
        p20.pick_up_tip()
        p20.aspirate(5, s.bottom(z=asp_height))
        p20.air_gap(airgap)
        ctx.delay(seconds=asp_delay)
        p20.dispense(5+airgap, d)
        p20.mix(mix_reps, 15, d)
        p20.blow_out()
        p20.drop_tip()
