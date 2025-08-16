import os, sys, pathlib

sys.path.append(str(pathlib.Path(__file__).resolve().parents[2] / 'packages' / 'electrofit' / 'src'))

import electrofit.core.process_initial_structure as core


def test_process_initial_structure_bcc_creates_gmx_files(tmp_path, monkeypatch):
    work = tmp_path / 'work'
    work.mkdir()
    mol2 = work / 'lig.mol2'
    mol2.write_text("@<TRIPOS>MOLECULE\nLIG\n 1 0 0 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n1 C1 0 0 0 C 1 LIG 0.0\n")
    equiv = work / 'equiv_groups.json'; equiv.write_text('{}')

    # Monkeypatch run_acpype symbol used inside core module to fabricate outputs
    def _fake_run_acpype(mol2_file, net_charge, scratch_dir, atom_type, charges="bcc"):
        base = os.path.splitext(os.path.basename(mol2_file))[0]
        ac_dir = os.path.join(scratch_dir, f"{base}.acpype")
        os.makedirs(ac_dir, exist_ok=True)
        for name, content in {
            f'{base}_GMX.gro': ';gro',
            f'{base}_GMX.itp': ';itp',
            f'{base}_GMX.top': '[ system ]\nX\n\n[ molecules ]\nX 1\n',
            f'posre_{base}.itp': '[ position_restraints ]\n',
        }.items():
            (pathlib.Path(ac_dir)/name).write_text(content)
    monkeypatch.setattr(core, 'run_acpype', _fake_run_acpype)

    old = os.getcwd(); os.chdir(work)
    try:
        core.process_initial_structure(
            molecule_name='lig',
            mol2_file=str(mol2),
            base_scratch_dir=str(tmp_path/'scratch_base'),
            additional_input=[str(equiv)],
            net_charge=0,
            residue_name='LIG',
            adjust_sym=False,
            ignore_sym=False,
            atom_type='gaff2',
            exit_screen=False,
            protocol='bcc'
        )
    finally:
        os.chdir(old)

    names = os.listdir(work)
    assert 'lig_GMX.gro' in names
    assert 'lig_GMX.itp' in names
    assert 'lig_GMX.top' in names
    assert any(n.startswith('posre_lig') and n.endswith('.itp') for n in names)
