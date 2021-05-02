import * as React from "react";
import { ICollection, IObjective, IParticleSelection } from "../../interfaces";
import { ButtonWithUpload } from "../ButtonWithUpload";
import { runMSO } from "../../utils/api";
import { FormWrapper } from "./FormWrapper";
import { UseStructureInputAddon } from "../UseStructureInputAddon";
import { Alert } from "react-bootstrap";
import { CurveEditorModal } from "../CurveEditorModal";

export const MSOForm = ({
  availableObjectives,
  addCollection,
  selection,
  loading,
  setLoading,
}: {
  availableObjectives?: IObjective[];
  addCollection(collection: ICollection): void;
  selection: IParticleSelection;
  loading: boolean;
  setLoading(loading: boolean): void;
}) => {
  const [name, setName] = React.useState<string>("MSO");
  const [msoStartingStructure, setMsoStartingStructure] = React.useState<string>("");
  const [nrOfParticles, setNrOfParticles] = React.useState<number>(50);
  const [nrOfIterations, setNrOfIterations] = React.useState<number>(10);
  const [nrOfSwarms, setNrOfSwarms] = React.useState<number>(5);
  const [vMin, setVMin] = React.useState<number>(-0.6);
  const [vMax, setVMax] = React.useState<number>(0.6);
  const [inertiaWeight, setInertiaWeight] = React.useState<number>(0.9);
  const [phi1, setPhi1] = React.useState<number>(2.0);
  const [phi2, setPhi2] = React.useState<number>(2.0);
  const [phi3, setPhi3] = React.useState<number>(2.0);
  const [selectedObjective, setSelectedObjective] = React.useState<any | null>();
  const [objectives, setObjectives] = React.useState<IObjective[]>([]);
  const [desirabilityCurveObjective, setDesirabilityCurveObjective] = React.useState<IObjective | null>(null);

  return (
    <FormWrapper
      title="Molecular Swarm Optimization"
      loading={loading}
      setLoading={setLoading}
      onSubmit={async () => {
        const data = await runMSO({
          structure: msoStartingStructure,
          iterations: nrOfIterations,
          num_swarms: nrOfSwarms,
          num_part: nrOfParticles,
          v_min: vMin,
          v_max: vMax,
          inertia_weight: inertiaWeight,
          phi1,
          phi2,
          phi3,
          objectives,
        });

        addCollection({ data, name });
      }}
    >
      {desirabilityCurveObjective ? (
        <CurveEditorModal
          open={desirabilityCurveObjective != null}
          initialPoints={desirabilityCurveObjective.desirability}
          setOpen={(open) => (open ? null : setDesirabilityCurveObjective(null))}
          onSave={(points) => {
            setObjectives(
              objectives.map((o) => (o === desirabilityCurveObjective ? { ...o, desirability: points } : o))
            );
            setDesirabilityCurveObjective(null);
          }}
        />
      ) : null}
      {availableObjectives ? (
        <>
          <div className="form-group">
            <label htmlFor="msoObjectivesSelect">Objectives</label>
            <div className="input-group input-group-sm">
              <select
                className="form-control"
                id="msoObjectivesSelect"
                value={selectedObjective?.name}
                onChange={(e) => {
                  setSelectedObjective(availableObjectives.find((o) => o.name === e.currentTarget.value));
                }}
              >
                <option value="">Choose...</option>
                {availableObjectives.map((o) => (
                  <option key={o.name} value={o.name}>{o.name}</option>
                ))}
              </select>
              <input
                type="number"
                className="form-control-plaintext form-control-sm"
                required
                value={objectives.length}
                onChange={() => null}
                min={1}
                style={{
                  flex: "0 0 1px",
                }}
              />
              <div className="input-group-append">
                <button
                  disabled={!selectedObjective}
                  className="btn btn-outline-secondary"
                  type="button"
                  onClick={() => setObjectives([...objectives, selectedObjective])}
                >
                  Add
                </button>
              </div>
            </div>
            <small id="msoObjectivesSelect" className="form-text text-muted">
              Objectives define the optimization goal of MSO
            </small>
          </div>
          <div className="form-group" style={{ overflowX: "hidden" }}>
            {objectives.map((objective) => (
              <details key={objective.name} className="mb-1" open={Boolean(objective.additional_args)}>
                <summary className="text-truncate">
                  {objective.name}{" "}
                  <div className="btn-group btn-group-sm ml-2 mr-2 float-right" role="group">
                    <button
                      type="button"
                      className="btn btn-light"
                      title="Adjust desirability curve"
                      onClick={() => {
                        setDesirabilityCurveObjective(objective);
                      }}
                    >
                      <i className="fas fa-fw fa-bezier-curve" />
                    </button>
                    <button
                      type="button"
                      className="btn btn-danger"
                      title="Delete objective"
                      onClick={() => {
                        setObjectives(objectives.filter((o) => o !== objective));
                      }}
                    >
                      <i className="fas fa-fw fa-times" />
                    </button>
                  </div>
                </summary>
                <small title={objective.description} className="text-truncate">
                  {objective.description}
                </small>
                <div className="d-flex" title="Relative weight of objective">
                  {/* <label for="formControlRange">Example Range input</label> */}
                  <input
                    type="range"
                    className="form-control-range mr-2"
                    min={1}
                    max={100}
                    step={1}
                    value={objective.weight}
                    onChange={(e) => {
                      setObjectives(
                        objectives.map((o) => (o === objective ? { ...o, weight: e.currentTarget.valueAsNumber } : o))
                      );
                    }}
                  />
                  {objective.weight}%
                </div>
                {Object.entries(objective.additional_args || {}).map(([key, value]) => {
                  const setValue = (value: string) =>
                    setObjectives(
                      objectives.map((o) =>
                        o === objective
                          ? {
                              ...o,
                              additional_args: {
                                ...o.additional_args,
                                [key]: value,
                              },
                            }
                          : o
                      )
                    );

                  return (
                    <div key={key} className="input-group input-group-sm">
                      {key === "query" ? (
                        <UseStructureInputAddon value={value} selection={selection} setValue={setValue} />
                      ) : null}
                      <input
                        key={key}
                        type="text"
                        className="form-control form-control-sm"
                        placeholder={key}
                        value={value}
                        required={true}
                        onChange={(e) => setValue(e.currentTarget.value)}
                      />
                    </div>
                  );
                })}
                <hr />
              </details>
            ))}
          </div>
          <div className="form-group">
            <label htmlFor="msoNameInput">Name</label>
            <input
              type="text"
              className="form-control form-control-sm"
              id="msoNameInput"
              required
              value={name}
              onChange={(e) => setName(e.currentTarget.value)}
            />
          </div>
          <div className="form-group">
            <label htmlFor="startingStructureInput">Starting structure</label>
            <div className="input-group input-group-sm">
              <UseStructureInputAddon
                value={msoStartingStructure}
                selection={selection}
                setValue={setMsoStartingStructure}
              />
              <input
                type="text"
                className="form-control form-control-sm"
                id="startingStructureInput"
                aria-describedby="startingStructureInputHelp"
                required
                value={msoStartingStructure}
                onChange={(e) => setMsoStartingStructure(e.currentTarget.value)}
              />
            </div>
            <small id="startingStructureInputHelp" className="form-text text-muted"></small>
          </div>
          <details>
            <summary>Advanced Settings</summary>
            <div className="form-group">
              <label htmlFor="nrOfSwarmsInput">Swarms</label>
              <input
                type="number"
                className="form-control form-control-sm"
                id="nrOfSwarmsInput"
                value={nrOfSwarms}
                required
                min={1}
                onChange={(e) => setNrOfSwarms(e.currentTarget.valueAsNumber)}
              />
            </div>
            <div className="form-group">
              <label htmlFor="nrOfParticlesInput">Particles</label>
              <input
                type="number"
                className="form-control form-control-sm"
                id="nrOfParticlesInput"
                value={nrOfParticles}
                required
                min={1}
                onChange={(e) => setNrOfParticles(e.currentTarget.valueAsNumber)}
              />
            </div>
            <div className="form-group">
              <label htmlFor="nrOfIterationsInput">Iterations</label>
              <input
                type="number"
                className="form-control form-control-sm"
                id="nrOfIterationsInput"
                value={nrOfIterations}
                required
                min={1}
                onChange={(e) => setNrOfIterations(e.currentTarget.valueAsNumber)}
              />
            </div>
            <div className="form-group">
              <label htmlFor="vMinInput">Min Velocity</label>
              <input
                type="number"
                className="form-control form-control-sm"
                id="vMinInput"
                value={vMin}
                required
                step={0.1}
                onChange={(e) => setVMin(e.currentTarget.valueAsNumber)}
              />
            </div>
            <div className="form-group">
              <label htmlFor="vMaxInput">Max Velocity</label>
              <input
                type="number"
                className="form-control form-control-sm"
                id="vMaxInput"
                value={vMax}
                required
                step={0.1}
                onChange={(e) => setVMax(e.currentTarget.valueAsNumber)}
              />
            </div>
            <div className="form-group">
              <label htmlFor="inertiaWeightInput">Inertia Weight</label>
              <input
                type="number"
                className="form-control form-control-sm"
                id="inertiaWeightInput"
                value={inertiaWeight}
                required
                step={0.1}
                onChange={(e) => setInertiaWeight(e.currentTarget.valueAsNumber)}
              />
            </div>
            <div className="form-group">
              <label htmlFor="phi1Input">Phi 1</label>
              <input
                type="number"
                className="form-control form-control-sm"
                id="phi1Input"
                value={phi1}
                required
                step={0.1}
                onChange={(e) => setPhi1(e.currentTarget.valueAsNumber)}
              />
            </div>
            <div className="form-group">
              <label htmlFor="phi2Input">Phi 2</label>
              <input
                type="number"
                className="form-control form-control-sm"
                id="phi2Input"
                value={phi2}
                required
                step={0.1}
                onChange={(e) => setPhi2(e.currentTarget.valueAsNumber)}
              />
            </div>
            <div className="form-group">
              <label htmlFor="phi3Input">Phi 3</label>
              <input
                type="number"
                className="form-control form-control-sm"
                id="phi3Input"
                value={phi3}
                required
                step={0.1}
                onChange={(e) => setPhi3(e.currentTarget.valueAsNumber)}
              />
            </div>
          </details>
          <div className="text-right">
            <ButtonWithUpload
              loading={loading}
              text="Run MSO"
              onUploadResult={(value) => {
                if (value) {
                  addCollection({ data: JSON.parse(value), name: "MSO" });
                  return true;
                }
              }}
            />
          </div>
        </>
      ) : (
        <Alert variant="info">
          <p>No objectives are available.</p>
        </Alert>
      )}
    </FormWrapper>
  );
};
