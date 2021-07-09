from opentrons import protocol_api
from opentrons.types import Point
import random

metadata = {
    "apiLevel": "2.8",
    "protocolName": "DNA_extraction_beads",
    "author": "REM",
    "description": "DNA extraction 96w plate",
}


def run(protocol: protocol_api.ProtocolContext):
    # columns
    column_no = 4
    used_cols = [f"A{i+1}" for i in range(column_no)]
    # Number of mixing times for lysate + beads
    lysate_mix_reps = 25
    # Number of mixing times for elution buffer
    eb_mix_reps = 25
    # lysate and bead volumes
    [v_sample, v_beads] = [50, 20]

    # ethanol volume in uL
    v_etoh = 120
    # final elution buffer
    v_EB = 50
    # binding time in minutes
    bead_binding_time = 10
    mag_binding_time = 5
    bead_to_magnet_time_in_s = 30
    air_dry_at_room_t = 5  # 5 minutes
    # beads: 1st row
    bead_tubes = protocol.load_labware(
        "opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap",
        7,
        label="bead-tubes",
    ).wells()[:2]
    # ethanol
    reagents = protocol.load_labware(
        "nest_12_reservoir_15ml", 2, label="etoh-eb trough"
    )
    etoh_loc1 = reagents["A2"]
    etoh_loc2 = reagents["A3"]

    # elution buffer
    eb_loc = reagents["A5"]
    # tip racks
    n_tiprack = 3

    # multi-channel tiprack
    slots = [3, 6, 9]
    assert n_tiprack <= len(slots)
    m_tipacks = [
        protocol.load_labware("opentrons_96_tiprack_300ul", i)
        for i in slots[:n_tiprack]
    ]

    # single-channel tiprack
    s_tiprack = protocol.load_labware("opentrons_96_tiprack_300ul", 8)

    # reaction plate (single one that will be moved by user)
    central_plate = protocol.load_labware(
        "96w_pcr_plate2", 5, label="central-plate"
    )
    mag_mod: protocol_api.MagneticModuleContext = protocol.load_module(
        "magnetic module", 1
    )
    rxn_plate = mag_mod.load_labware("96w_pcr_plate2", label="magnet-plate")
    elution_plate = protocol.load_labware(
        "96w_pcr_plate2", 4, label="elution-plate"
    )

    # pipettes
    # multi-channel
    multi_300: protocol_api.InstrumentContext = protocol.load_instrument(
        "p300_multi_gen2", mount="right", tip_racks=m_tipacks
    )
    tip300_count = 0
    tip300_max = 12 * n_tiprack

    def multi_tip_pick_up():
        """Refillable tip pickup.
        If no tips left, protocol will
        pause and ask for refill."""
        nonlocal tip300_count
        if tip300_count == tip300_max:
            protocol.pause(
                f"Refill 300ul tipracks in slots {slots} before resuming."
            )
            multi_300.reset_tipracks()
            tip300_count = 0
        tip300_count += 1
        multi_300.pick_up_tip()

    # single channel
    single_300: protocol_api.InstrumentContext = protocol.load_instrument(
        "p300_single_gen2", mount="left", tip_racks=[s_tiprack]
    )

    # lights on
    protocol.set_rail_lights(False)
    protocol.set_rail_lights(True)

    # step 1 magnet
    mag_mod.disengage()
    protocol.pause("Put plate on slot 5 (not magnet) and start.")
    # step 2 add beads
    protocol.comment("Bead transfer")
    single_300.flow_rate.aspirate = 120
    single_300.flow_rate.dispense = 120

    for wells in central_plate.columns()[:column_no]:
        single_300.distribute(
            v_beads,
            bead_tubes,
            wells,
            air_gap=20,
            new_tip="always",
            mix_before=(5, 200),
            blow_out=True,
            blow_out_location="destination_well",
        )

    protocol.comment("Lysate mixing")
    # step 3 mix with multi-channel
    multi_300.flow_rate.aspirate = 67
    multi_300.flow_rate.dispense = 120

    # 2 ok for 50uL+20uL, 4 for 120uL total
    multi_300.well_bottom_clearance.aspirate = 3
    multi_300.well_bottom_clearance.dispense = 2

    # mix and return tips
    for col in used_cols:
        multi_300.pick_up_tip()
        multi_300.mix(
            lysate_mix_reps,
            min((v_beads + v_sample) * 0.7, multi_300.max_volume),
            central_plate[col],
        )
        multi_300.blow_out(central_plate[col].top())
        multi_300.return_tip()
    multi_300.reset_tipracks()

    # step 4 bead binding
    protocol.delay(
        minutes=bead_binding_time,
        msg=f"RT: {str(bead_binding_time)} minutes for\
             DNA binding to beads.",
    )
    protocol.pause(
        "Place sample plate on magnetic module.\
        Place a trash plate on slot 5"
    )

    # step 5 engage magnet
    mag_mod.engage(height=15)
    protocol.delay(
        minutes=mag_binding_time,
        msg=f"Wait {str(mag_binding_time)} minutes for\
             beads to bind to magnet.Make sure beads\
                  are visibly sticking to the wall!",
    )

    # step 6 transfer EDM supernatant to waste
    def mag_x_location(col_no: int) -> int:
        """Return x coordinate of magnet
        as a function of column no that can be [0,11]"""
        # left if even and right if odd
        return -1.1 if (col_no + 1) % 2 == 0 else 1.1

    def transfer_supernatant(vol, src_pl, dest_pl):
        """Will transfer supernatant avoiding to take bead pellets.
        Picks up tips, transfers with air gap of 20uL
        and then drops tips"""
        for i, col in enumerate(used_cols):
            multi_tip_pick_up()
            # offset x location dependening on magnet location
            multi_300.transfer(
                vol,
                src_pl[col].bottom(z=1).move(Point(x=-mag_x_location(i))),
                dest_pl[col],
                air_gap=20,
                new_tip="never",
            )
            multi_300.drop_tip()

    v_asp_sn = (v_sample + v_beads) * 1.4
    multi_300.flow_rate.aspirate = 30
    multi_300.flow_rate.dispense = 300
    # multi_300.well_bottom_clearance.dispense = 5
    # multi_300.well_bottom_clearance.aspirate = 1.7
    transfer_supernatant(v_asp_sn, rxn_plate, central_plate)
    # step 7 disengage magnet
    mag_mod.disengage()
    multi_300.flow_rate.aspirate = 150
    multi_300.flow_rate.dispense = 300
    multi_300.well_bottom_clearance.dispense = 6

    # step 8: ethanol wash. step 9: mix
    def mix_at_custom_speed_at_random_offsets(
        reps: int,  # repetitions
        vol: float,  # volumes
        dest: Point,  # destination point
        asp_r: float,  # aspirationrate
        disp_r: float = None,  # dispense rate defaults to equal asp_r
        x_limits=(-1.1, 1.1),  # min-max pipette locations
        y_limits=(-0.8, 0.8),
        z_limits=(0.2, 3),
    ):
        """Wrapper of mixing functionality.
        We customize aspiration and dispense speeds.
        reps: repetitions
        vol: volume of mixing
        dest: destination of mixing
        asp_r,disp_r: multiplication factor
        of base flow rate.
        If 4, then flow rate will be 4 times base flow rate
        x_limits, y_limits, z_limits: range of locations
        for pipette within well"""
        if not disp_r:
            disp_r = asp_r
        for i in range(reps):
            x = round(random.uniform(x_limits[0], x_limits[1]), 1)
            y = round(random.uniform(y_limits[0], y_limits[1]), 1)
            z = round(random.uniform(z_limits[0], z_limits[1]), 1)
            offset = Point(x, y, z)
            multi_300.aspirate(vol, dest.move(offset), asp_r)
            multi_300.dispense(vol, dest.move(offset), disp_r)

    tot_v_etoh = 0
    for wash_cycle in range(2):
        for i, col in enumerate(used_cols):
            tot_v_etoh += v_etoh
            # transfer from second column/third column
            etoh_source = etoh_loc1 if tot_v_etoh // 1000 < 1 else etoh_loc2
            multi_300.pick_up_tip()
            # transfer ethanol to reaction plate
            multi_300.transfer(
                v_etoh,
                etoh_source,
                rxn_plate[col],
                rate=0.5,
                air_gap=30,
                new_tip="never",
            )
            base_h = 1  # min height from base in mm
            # mix at bead pellet
            mix_at_custom_speed_at_random_offsets(
                5,  # repetitions
                v_etoh,  # volume
                rxn_plate[col].bottom(base_h),  # destination
                5,  # times base flow rate
                x_limits=(
                    mag_x_location(i),
                    mag_x_location(i),
                ),  # focus on bead pellet
            )
            # generalized mix
            mix_at_custom_speed_at_random_offsets(
                10, v_etoh, rxn_plate[col].bottom(base_h), 4
            )
            protocol.delay(seconds=2)
            multi_300.blow_out(rxn_plate[col].top(z=-1))
            multi_300.return_tip()
        # step 10 engage magnet
        mag_mod.engage(height=15)
        multi_300.reset_tipracks()

        # step 11
        protocol.delay(
            seconds=bead_to_magnet_time_in_s,
            msg=f"Turn on heatblock at 60°C.\
                 Wait {bead_to_magnet_time_in_s} until magnet acts properly.\
                 make sure it is sticking to wall, or pause manually.\
                      If needed remove tips from trash.",
        )
        # step 12
        transfer_supernatant(
            min(1.2 * v_etoh, multi_300.max_volume - 30),
            rxn_plate,
            central_plate,
        )
        if wash_cycle == 0:
            # step 13 disengage magnet
            mag_mod.disengage()
    # step 14 air dry and dry at heating block
    protocol.delay(
        minutes=air_dry_at_room_t,
        msg=f"Air dry {air_dry_at_room_t} minutes.",
    )
    protocol.pause(
        "Put plate on thermoblock for 2 min to finish drying the wells.\
        Continue when ethanol has dried up or been removed.\
            Put a clean plate on slot 5 to collect samples."
    )
    # step 15 disengage mag mod
    mag_mod.disengage()
    # step 16 EB mixing
    for col in used_cols:
        multi_tip_pick_up()
        multi_300.transfer(
            v_EB, eb_loc, rxn_plate[col], rate=0.5, new_tip="never", air_gap=20
        )
        mix_at_custom_speed_at_random_offsets(
            eb_mix_reps, v_EB * 0.7, rxn_plate[col].bottom(1), 4
        )
        multi_300.drop_tip()
    # step 17 increase yield by transferring to heating block
    protocol.pause(
        "Put plate 5min at 60°C with 300rpm shaking\
             and put back onto the magnet before continuing.\
                  Make sure to have a new plate on slot\
                       5 for final DNA samples"
    )
    # step 18 engage mag mod to pull pure dna
    mag_mod.engage(height=14)
    protocol.delay(
        seconds=bead_to_magnet_time_in_s,
        msg=f"{bead_to_magnet_time_in_s} seconds for magnet action on beads",
    )
    multi_300.flow_rate.aspirate = 20
    multi_300.flow_rate.dispense = 150
    # step 19
    transfer_supernatant(v_EB * 1.2, rxn_plate, elution_plate)
