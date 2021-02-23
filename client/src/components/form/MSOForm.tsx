import * as React from "react";
import { ICollection, IParticleSelection } from "../../interfaces";
import { ButtonWithUpload } from "../ButtonWithUpload";
import { runMSO } from "../../utils/api";
import { FormWrapper } from "./FormWrapper";
import { StructureImage } from "../StructureImage";
import { UseStructureInputAddon } from "../UseStructureInputAddon";

export const MSOForm = ({
  addCollection,
  selection,
  loading,
  setLoading,
}: {
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

  return (
    <FormWrapper
      title="Run MSO from starting structure"
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
        });

        addCollection({ data, name });
      }}
    >
      <div className="form-group">
        <label htmlFor="msoNameInput">Name</label>
        <input
          type="text"
          className="form-control form-control-sm"
          id="msoNameInput"
          value={name}
          onChange={(e) => setName(e.currentTarget.value)}
        />
      </div>
      <div className="form-group">
        <label htmlFor="nrOfSwarmsInput">Swarms</label>
        <input
          type="number"
          className="form-control form-control-sm"
          id="nrOfSwarmsInput"
          value={nrOfSwarms}
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
          min={1}
          onChange={(e) => setNrOfIterations(e.currentTarget.valueAsNumber)}
        />
      </div>
      {/* 
  v_min?: number;
  v_max?: number;
  inertia_weight?: number;
  phi1?: number;
  phi2?: number;
  phi3?: number; */}
      <div className="form-group">
        <label htmlFor="vMinInput">Min Velocity</label>
        <input
          type="number"
          className="form-control form-control-sm"
          id="vMinInput"
          value={vMin}
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
          step={0.1}
          onChange={(e) => setPhi3(e.currentTarget.valueAsNumber)}
        />
      </div>
      <div className="form-group">
        <label htmlFor="startingStructureInput">Starting structure</label>
        <div style={{ height: 100 }}>
          {msoStartingStructure ? (
            <StructureImage
              style={{ height: "100%", margin: "0 auto", display: "block" }}
              structure={msoStartingStructure}
            />
          ) : null}
        </div>
        <div className="input-group input-group-sm">
          <UseStructureInputAddon selection={selection} setValue={setMsoStartingStructure} />
          <input
            type="text"
            className="form-control form-control-sm"
            id="startingStructureInput"
            aria-describedby="startingStructureInputHelp"
            value={msoStartingStructure}
            onChange={(e) => setMsoStartingStructure(e.currentTarget.value)}
          />
        </div>
        <small id="startingStructureInputHelp" className="form-text text-muted"></small>
      </div>
      <div className="text-right">
        <ButtonWithUpload
          loading={loading}
          disabled={!msoStartingStructure}
          text="Run MSO"
          onUploadResult={(value) => {
            if (value) {
              addCollection({ data: JSON.parse(value), name: "MSO" });
              return true;
            }
          }}
        />
      </div>
    </FormWrapper>
  );
};
