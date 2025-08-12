import typer
from electrofit.workflows import (
    step0_setup, step1_initial_processing, step2_setup_sim_dir,
    step3_start_sim, step4_extract_conforms, step5_process_conforms,
    step6_extract_average_charges, step7_setup_final_sim, step8_start_final_sim,
)

app = typer.Typer(help="Electrofit core CLI")

@app.command("run")
def run_all():
    step0_setup.run()
    step1_initial_processing.run()
    step2_setup_sim_dir.run()
    step3_start_sim.run()
    step4_extract_conforms.run()
    step5_process_conforms.run()
    step6_extract_average_charges.run()
    step7_setup_final_sim.run()
    step8_start_final_sim.run()

def main():
    app()