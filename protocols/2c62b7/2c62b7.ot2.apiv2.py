from opentrons import types
from opentrons import protocol_api
import math

metadata = {
    'protocolName': 'Omega Bio-tek Mag-Bind Environmental DNA 96 Kit',
    'author': 'Sakib <sakib.hossain@opentrons.com>',
    'description': 'Custom Protocol Request',
    'apiLevel': '2.9'
}


def run(ctx):

    [debug, samples, m300_mount, m20_mount, tip_type, mm2_vol, vhb_vol,
        elution_buffer_vol, settling_time] = get_values(  # noqa: F821
        "debug", "samples", "m300_mount", "m20_mount", "tip_type", "mm2_vol",
        "vhb_vol", "elution_buffer_vol", "settling_time")

    cols = math.ceil(samples/8)

    tiprack_type = {
        'standard': 'opentrons_96_tiprack_300ul',
        'filter': 'opentrons_96_filtertiprack_200ul'
        }

    # Load Labware/Modules
    temp_mod = ctx.load_module('temperature module gen2', 3)
    temp_plate = temp_mod.load_labware(
                    'opentrons_96_aluminumblock_nest_wellplate_100ul')
    mag_mod = ctx.load_module('magnetic module gen2', 1)
    mag_plate = mag_mod.load_labware('nest_96_wellplate_2ml_deep')
    tipracks = [ctx.load_labware(tiprack_type[tip_type], slot) for slot in
                [4, 7, 8, 9, 10]]
    tip_isolator = ctx.load_labware(tiprack_type[tip_type], 11, 'Tip Isolator')
    res1 = ctx.load_labware('nest_12_reservoir_15ml', 5)
    res2 = ctx.load_labware('nest_12_reservoir_15ml', 2)
    dna_plate = ctx.load_labware('nest_96_wellplate_100ul_pcr_full_skirt', 6)
    trash = ctx.loaded_labwares[12]['A1']

    # Load Pipettes
    m300 = ctx.load_instrument('p300_multi_gen2', m300_mount,
                               tip_racks=tipracks)
    max_tip_volume = tipracks[0]['A1'].geometry.max_volume

    # Reagents
    # Splitting columns for an even 12 column transfer
    # based on volume total
    mm2 = [well for well in res1.wells()[:3] for i in range(4)]
    vhb = [well for well in res1.wells()[4:8] for i in range(3)]
    etoh1 = [well for well in res2.wells()[:4] for i in range(3)]
    etoh2 = [well for well in res2.wells()[4:8] for i in range(3)]
    elution_buffer = res1.wells()[11]
    elution_wells = temp_plate.rows()[0][:cols]

    # Helper Functions
    def debug_mode(msg, debug_setting=debug):
        if debug_setting == "True":
            ctx.pause(msg)

    def supernatant_removal(vol, src, dest, side=-1):
        m300.flow_rate.aspirate = 20
        while vol >= max_tip_volume:
            m300.aspirate(
                max_tip_volume, src.bottom().move(
                    types.Point(x=side, y=0, z=0.5)))
            m300.dispense(max_tip_volume, dest)
            vol -= max_tip_volume

        if vol < max_tip_volume:
            m300.aspirate(vol, src.bottom().move(
                        types.Point(x=side, y=0, z=0.5)))
            m300.dispense(vol, dest)
        m300.flow_rate.aspirate = 50

    def reset_flow_rates():
        m300.flow_rate.aspirate = 94
        m300.flow_rate.dispense = 94

    def pick_up(pip, loc=None):
        """Function that can be used instead of .pick_up_tip() that will pause
        robot when robot runs out of tips, prompting user to replace tips
        before resuming"""
        try:
            if loc:
                pip.pick_up_tip(loc)
            else:
                pip.pick_up_tip()
        except protocol_api.labware.OutOfTipsError:
            pip.home()
            ctx.pause("Replace the tips")
            pip.reset_tipracks()
            pip.pick_up_tip()

    def tip_mix(well, vol, reps, park_tip=False, tip_loc=None, tip_map=None,
                asp_speed=94, disp_speed=94):

        m300.flow_rate.aspirate = asp_speed
        m300.flow_rate.dispense = disp_speed

        if not m300.has_tip:
            if tip_loc:
                pick_up(m300, tip_loc)
            else:
                pick_up(m300)
        ctx.comment('Mixing from the middle')
        m300.mix(reps, vol, well.bottom(z=4))
        ctx.comment('Mixing from the bottom')
        m300.mix(reps, vol, well.bottom())
        ctx.comment('Mixing from the middle')
        m300.mix(reps, vol, well.bottom(z=4))
        m300.blow_out()
        m300.touch_tip()

        if not park_tip:
            m300.drop_tip()
        elif park_tip:
            m300.drop_tip(tip_isolator.columns()[mag_plate_wells[well]][0])
        reset_flow_rates()

    def bind(delay=settling_time):
        if mag_mod.status != 'engaged':
            mag_mod.engage()
        ctx.delay(minutes=delay, msg=f'''Incubating on MagDeck for
                  {delay} minute(s).''')

    def remove(vol, src, dest=trash, use_park_tip=True):
        if use_park_tip:
            pick_up(m300, tip_isolator.columns()[mag_plate_wells[src]][0])
        elif not use_park_tip:
            pick_up(m300)
        supernatant_removal(vol=vol, src=src, dest=dest)
        m300.drop_tip()

    def wash(vol, src, dest):
        m300.flow_rate.dispense = 200
        pick_up(m300)
        m300.transfer(vol, src, dest, new_tip='never')
        m300.drop_tip()
        m300.flow_rate.dispense = 94

    def etoh_wash(reservoir, park=True):
        # Steps 24-25
        # 70% Ethanol Wash
        debug_mode(msg="Debug: Wash with 70% Ethanol")
        for src, dest in zip(reservoir, mag_plate_wells):
            wash(500, src, dest)
        # Tip Mix (Vortex)
        debug_mode(msg="Debug: Tip Mixing (Vortex)")
        for well in mag_plate_wells:
            tip_mix(well, 300, 10, park_tip=park,
                    tip_loc=tip_isolator.columns()[mag_plate_wells[well]][0],
                    tip_map=mag_plate_wells, asp_speed=188, disp_speed=188)

        # Steps 26-27
        debug_mode(msg=f'''Debug: Engage Magnet for {settling_time} minutes
                        and then remove supernatant''')
        bind()
        for well in mag_plate_wells:
            remove(500, well, use_park_tip=False)
        if mag_mod.status == 'engaged':
            mag_mod.disengage()

    # Volume Tracking
    class VolTracker:
        def __init__(self, labware, well_vol, pip_type='single',
                     mode='reagent', start=0, end=12):
            self.labware_wells = dict.fromkeys(labware.wells()[start:end], 0)
            self.well_vol = well_vol
            self.pip_type = pip_type
            self.mode = mode
            self.start = start
            self.end = end

        def tracker(self, vol):
            '''tracker() will track how much liquid
            was used up per well. If the volume of
            a given well is greater than self.well_vol
            it will remove it from the dictionary and iterate
            to the next well which will act as the reservoir.'''
            well = next(iter(self.labware_wells))
            if self.labware_wells[well] >= self.well_vol:
                del self.labware_wells[well]
                well = next(iter(self.labware_wells))
            if self.pip_type == 'multi':
                self.labware_wells[well] = self.labware_wells[well] + vol*8
            elif self.pip_type == 'single':
                self.labware_wells[well] = self.labware_wells[well] + vol
            if self.mode == 'waste':
                ctx.comment(f'''{well}: {int(self.labware_wells[well])} uL of
                            total waste''')
            else:
                ctx.comment(f'''{int(self.labware_wells[well])} uL of liquid
                            used from {well}''')
            return well

    # Track Reagent Volumes
    water = VolTracker(res2, 14400, 'multi', 'reagent', start=9, end=12)

    # Wells
    # mag_plate_wells = mag_plate.rows()[0]
    mag_plate_wells = {well: column for well, column in zip(
                       mag_plate.rows()[0][:cols], range(cols))}
    dna_plate_wells = dna_plate.rows()[0][:cols]

    # Protocol Steps

    # Add 300 uL of Water to each well in tip isolator
    ctx.comment('''Transferring 300 uL of water to
                each well in the tip isolator''')
    pick_up(m300, tip_isolator['A1'])
    for col in tip_isolator.rows()[0]:
        m300.transfer(300, water.tracker(300), col, new_tip='never')

    # Step 14
    # Transfer Master Mix 2 (XP1 Buffer + Mag-Bind® Particles RQ) to Mag Plate
    debug_mode(msg="Debug: Transfer Master Mix 2 to Mag Plate (Step 14)")
    transfer_count = 0
    for mm, dest in zip(mm2, mag_plate_wells):
        if transfer_count == 3:
            transfer_count = 0
        if transfer_count == 0:
            if not m300.has_tip:
                pick_up(m300)
            m300.mix(10, 300, mm.bottom(z=4))
            if cols > 6:
                m300.mix(10, 300, mm2[6].bottom())
            m300.drop_tip()
        pick_up(m300)
        m300.aspirate(mm2_vol, mm)
        m300.dispense(mm2_vol, dest)
        m300.drop_tip()
        transfer_count += 1

    # Step 15
    # Incubate at Room Temp for 10 Minutes while mixing
    debug_mode(msg='''Debug: Incubate at Room Temperature while
               Mixing (Step 15)''')
    for well in mag_plate_wells:
        tip_mix(well, mm2_vol/2, 10, park_tip=True, tip_map=mag_plate_wells,
                asp_speed=94, disp_speed=94)

    for well in mag_plate_wells:
        tip_mix(well, mm2_vol/2, 5, park_tip=True,
                tip_loc=tip_isolator.columns()[mag_plate_wells[well]][0],
                tip_map=mag_plate_wells, asp_speed=94, disp_speed=94)

    # Steps 16-18
    debug_mode(msg=f'''Debug: Engage Magnet for {settling_time} minutes
                    and then remove supernatant (Steps 16-18)''')
    bind()
    for well in mag_plate_wells:
        remove(600, well, use_park_tip=True)
    if mag_mod.status == 'engaged':
        mag_mod.disengage()

    # Steps 19-20
    # Add VHB Buffer
    debug_mode(msg="Debug: Wash with VHB Buffer (Step 19)")
    for src, dest in zip(vhb, mag_plate_wells):
        wash(500, src, dest)
    # Tip Mix (Vortex)
    debug_mode(msg="Debug: Tip Mixing (Vortex) (Step 20)")
    for well in mag_plate_wells:
        tip_mix(well, 300, 10, park_tip=True, tip_map=mag_plate_wells,
                asp_speed=188, disp_speed=188)

    # Steps 21-23
    debug_mode(msg=f'''Debug: Engage Magnet for {settling_time} minutes and
               then remove supernatant (Steps 21-23)''')
    bind()
    for well in mag_plate_wells:
        remove(400, well, use_park_tip=False)
    if mag_mod.status == 'engaged':
        mag_mod.disengage()

    # Steps 24-28
    etoh_wash(etoh1)
    etoh_wash(etoh2, park=False)

    # Step 29
    debug_mode(msg='''Debug: Engaging Magnet for 1 minute and removing any
               supernatant (Step 29)''')
    bind(delay=1)
    for well in mag_plate_wells:
        remove(200, well, dest=trash, use_park_tip=False)

    # Step 30
    debug_mode(msg='''Debug: Engaging Magnet for 10 minutes to allow beads to
               dry (Step 30)''')
    bind(delay=10)

    # Step 31
    debug_mode(msg='''Debug: Transferring Elution Buffer to plate on
               temperature module (Step 31)''')
    m300.transfer(elution_buffer_vol, elution_buffer, elution_wells)
    debug_mode(msg='''Debug: Heating elution buffer on temperature module to
               70C (Step 31)''')
    temp_mod.set_temperature(70)
    debug_mode(msg='''Debug: Transferring elution buffer to sample wells
               (Step 31)''')
    for src, dest in zip(elution_wells, mag_plate_wells):
        m300.transfer(elution_buffer_vol, src, dest)

    # Step 32
    debug_mode(msg='''Debug: Tip Mixing (Vortex) (Step 32)''')
    for well in mag_plate_wells:
        tip_mix(well, elution_buffer_vol/2, 10, park_tip=False, tip_loc=None,
                tip_map=None, asp_speed=94, disp_speed=94)

    # Step 33
    debug_mode(msg='''Debug: Engaging Magnetic Module for 2 minutes to allow
               beads to settle (Step 33)''')
    bind(delay=2)

    # Step 34
    debug_mode(msg='''Debug: Transfer clear supernatant containing purified
               DNA to NEST 0.1 mL 96 Well PCR Plate (Step 34)''')
    for src, dest in zip(elution_wells, dna_plate_wells):
        pick_up(m300)
        m300.aspirate(elution_buffer_vol, src)
        m300.dispense(elution_buffer_vol, dest.bottom(z=5))
        m300.drop_tip()
