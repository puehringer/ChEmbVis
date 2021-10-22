import * as React from "react";
import { ICollection, IParticleSelection } from "../../interfaces";
import { ButtonWithUpload } from "../ButtonWithUpload";
import { getMatchedMolecularPairs, embedStructures, getStonedSelfies } from "../../utils/api";
import { FormWrapper } from "./FormWrapper";
import { UseStructureInputAddon } from "../UseStructureInputAddon";
import { useNameInput } from "../../utils/hooks";

export const StonedSelfiesForm = ({
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
  const [structure, setStructure] = React.useState<string>("");
  const [name, setName, nameInput] = useNameInput("stonedSelfiesNameInput", `SELFIES ${structure}`);
  const [substructure, setSubstructure] = React.useState<string>("");
  const [randomSamples, setRandomSamples] = React.useState<number>(1000);
  const [maxMutations, setMaxMutations] = React.useState<number>(5);

  return (
    <FormWrapper
      title={
        <>
          Stoned Selfies &nbsp;
          <a
            style={{ fontSize: "smaller" }}
            href="https://github.com/aspuru-guzik-group/stoned-selfies"
            target="_blank"
            rel="noreferrer"
          >
            <i className="fas fa-fw fa-info-circle" />
          </a>
        </>}
      loading={loading}
      setLoading={setLoading}
      onSubmit={async () => {
        const { structures } = await getStonedSelfies({
          structure,
          substructure,
          random_samples: randomSamples,
          max_mutations: maxMutations,
        });
        if (structures.length === 0) {
          throw Error("No stoned selfies found");
        }
        const embeddedCollection = await embedStructures({ structures, include_embedding: true });
        addCollection({
          name,
          ...embeddedCollection,
        });
        setName("");
      }}
    >
      {nameInput}
      <div className="mb-3">
        <label htmlFor="stonedSelfiesStructureInput">Structure</label>
        <div className="input-group input-group-sm">
          <UseStructureInputAddon value={structure} selection={selection} setValue={setStructure} />
          <input
            type="text"
            className="form-control form-control-sm"
            id="stonedSelfiesStructureInput"
            aria-describedby="stonedSelfiesStructureInput"
            value={structure}
            onChange={(e) => setStructure(e.currentTarget.value)}
          />
        </div>
      </div>
      <details>
        <summary>Advanced Settings</summary>
        <div className="mb-3">
          <label htmlFor="stonedSelfiesSubstructureInput">Substructure</label>
          <div className="input-group input-group-sm">
            <UseStructureInputAddon value={substructure} selection={selection} setValue={setSubstructure} />
            <input
              type="text"
              className="form-control form-control-sm"
              id="stonedSelfiesSubstructureInput"
              aria-describedby="stonedSelfiesSubstructureInput"
              value={substructure}
              onChange={(e) => setSubstructure(e.currentTarget.value)}
            />
          </div>
        </div>
        <div className="mb-3">
          <label htmlFor="stonedSelfiesRandomSamplesInput" className="form-label">Nr. of Random Samples: {randomSamples}</label>
          <input
            type="range"
            className="form-range"
            id="stonedSelfiesRandomSamplesInput"
            value={randomSamples}
            onChange={(e) => setRandomSamples(e.currentTarget.valueAsNumber)}
            min={500}
            max={50000}
          />
        </div>
        <div className="mb-3">
          <label htmlFor="stonedSelfiesMaxPermutationsInput" className="form-label">Maximum Permutations: {maxMutations}</label>
          <input
            type="range"
            className="form-range"
            id="stonedSelfiesMaxPermutationsInput"
            value={maxMutations}
            onChange={(e) => setMaxMutations(e.currentTarget.valueAsNumber)}
            min={1}
            max={20}
          />
        </div>
      </details>
      <div className="text-end">
        <ButtonWithUpload loading={loading} disabled={!structure} text="Compute" />
      </div>
    </FormWrapper>
  );
};
