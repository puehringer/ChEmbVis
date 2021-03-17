import * as React from "react";
import { ICollection, IParticleSelection } from "../../interfaces";
import { ButtonWithUpload } from "../ButtonWithUpload";
import { getMatchedMolecularPairs, embedStructures } from "../../utils/api";
import { FormWrapper } from "./FormWrapper";
import { UseStructureInputAddon } from "../UseStructureInputAddon";

export const MatchedMolecularPairsForm = ({
  addCollection,
  selection,
  loading,
  setLoading,
}: {
  addCollection(collection: ICollection): void;
  selection?: IParticleSelection;
  loading: boolean;
  setLoading(loading: boolean): void;
}) => {
  const [name, setName] = React.useState<string>("");
  const [structure, setStructure] = React.useState<string>("");
  const [substructure, setSubstructure] = React.useState<string>("");
  const [minVariableSize, setMinVariableSize] = React.useState<number>(0);
  const [maxVariableSize, setMaxVariableSize] = React.useState<number>(1);
  const [minConstantSize, setMinConstantSize] = React.useState<number>(0);
  const [minRadius, setMinRadius] = React.useState<number>(0);
  const [minPairs, setMinPairs] = React.useState<number>(0);

  return (
    <FormWrapper
      title={
        <>
          Matched Molecular Pairs &nbsp;
          <a
            style={{ fontSize: "smaller" }}
            href="https://github.com/rdkit/mmpdb#4-identify-possible-transforms"
            target="_blank"
            rel="noreferrer"
          >
            <i className="fas fa-fw fa-info-circle" />
          </a>
        </>
      }
      loading={loading}
      setLoading={setLoading}
      onSubmit={async () => {
        const {structures} = await getMatchedMolecularPairs({
          structure,
          substructure: substructure,
          min_variable_size: minVariableSize,
          max_variable_size: maxVariableSize,
          min_constant_size: minConstantSize,
          min_radius: minRadius,
          min_pairs: minPairs,
        });
        if(structures.length === 0) {
          throw Error('No matching pairs found');
        }
        const particles = await embedStructures({structures});
        addCollection({
          name: name || `MMP ${structure}`,
          data: particles,
        });
      }}
    >
      <div className="form-group">
        <label htmlFor="mmpStructureInput">Structure</label>
        <div className="input-group input-group-sm">
          <UseStructureInputAddon value={structure} selection={selection} setValue={setStructure} />
          <input
            type="text"
            className="form-control form-control-sm"
            id="mmpStructureInput"
            aria-describedby="mmpStructureInput"
            value={structure}
            onChange={(e) => setStructure(e.currentTarget.value)}
          />
        </div>
      </div>
      <details>
        <summary>Advanced Settings</summary>
        <div className="form-group">
          <label htmlFor="mmpNameInput">Name</label>
          <input
            type="text"
            className="form-control form-control-sm"
            id="mmpNameInput"
            placeholder={`MMP ${structure}`}
            value={name}
            onChange={(e) => setName(e.currentTarget.value)}
          />
        </div>
        <div className="form-group">
          <label htmlFor="minVariableSizeInput">Min Variable Size</label>
          <input
            type="number"
            className="form-control form-control-sm"
            id="minVariableSizeInput"
            value={minVariableSize}
            required
            min={0}
            onChange={(e) => setMinVariableSize(e.currentTarget.valueAsNumber)}
          />
        </div>
        <div className="form-group">
          <label htmlFor="maxVariableSizeInput">Max Variable Size</label>
          <input
            type="number"
            className="form-control form-control-sm"
            id="maxVariableSizeInput"
            value={maxVariableSize}
            required
            min={0}
            onChange={(e) => setMaxVariableSize(e.currentTarget.valueAsNumber)}
          />
        </div>
        <div className="form-group">
          <label htmlFor="minConstantSizeInput">Min Constant Size</label>
          <input
            type="number"
            className="form-control form-control-sm"
            id="minConstantSizeInput"
            value={minConstantSize}
            required
            min={0}
            onChange={(e) => setMinConstantSize(e.currentTarget.valueAsNumber)}
          />
        </div>
        <div className="form-group">
          <label htmlFor="minRadiusInput">Min Radius</label>
          <input
            type="number"
            className="form-control form-control-sm"
            id="minRadiusInput"
            value={minRadius}
            min={0}
            required
            onChange={(e) => setMinRadius(e.currentTarget.valueAsNumber)}
          />
        </div>
        <div className="form-group">
          <label htmlFor="minPairsInput">Min Pairs</label>
          <input
            type="number"
            className="form-control form-control-sm"
            id="minPairsInput"
            value={minPairs}
            min={0}
            required
            onChange={(e) => setMinPairs(e.currentTarget.valueAsNumber)}
          />
        </div>
        <div className="form-group">
          <label htmlFor="mmpSmartsStructureInput">SMARTS Substructure</label>
          <div className="input-group input-group-sm">
            <UseStructureInputAddon value={substructure} selection={selection} setValue={setSubstructure} />
            <input
              type="text"
              className="form-control form-control-sm"
              id="mmpSmartsStructureInput"
              aria-describedby="mmpSmartsStructureInput"
              value={substructure}
              onChange={(e) => setSubstructure(e.currentTarget.value)}
            />
          </div>
        </div>
      </details>
      <div className="text-right">
        <ButtonWithUpload loading={loading} disabled={!structure} text="Compute MMP" />
      </div>
    </FormWrapper>
  );
};
