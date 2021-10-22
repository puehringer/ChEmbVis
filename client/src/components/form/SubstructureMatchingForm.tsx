import * as React from "react";
import { ICollection, IParticleSelection } from "../../interfaces";
import { ButtonWithUpload } from "../ButtonWithUpload";
import { hasSubstructureMatch } from "../../utils/api";
import { FormWrapper } from "./FormWrapper";
import { UseStructureInputAddon } from "../UseStructureInputAddon";

export const SubstructureMatchingForm = ({
  collections,
  setCollections,
  selection,
  loading,
  setLoading,
}: {
  collections: ICollection[];
  setCollections(collections: ICollection[]): void;
  selection: IParticleSelection;
  loading: boolean;
  setLoading(loading: boolean): void;
}) => {
  const [smarts, setSmarts] = React.useState<string>("");

  const computeSubstructures = async () => {
    return Promise.all(
      collections.map(({ data }) =>
        hasSubstructureMatch(
          data.map(({ structure }) => structure),
          smarts
        )
      )
    ).then((results) => {
      setCollections(
        collections.map((c, i) => {
          return {
            ...c,
            data: c.data.map((p) => ({
              ...p,
              properties: {
                ...(p.properties || {}),
                [`Has ${smarts}`]: Boolean(results[i].validity[p.structure]),
                [`Count ${smarts}`]: results[i].counts[p.structure] ?? 0,
              },
            })),
            selection: null,
          };
        })
      );
    });
  };

  return (
    <FormWrapper
      title="Substructure Matching"
      loading={loading}
      setLoading={setLoading}
      onSubmit={computeSubstructures}
    >
      <div className="mb-3">
        <label htmlFor="smartsStructureInput">Substructure</label>
        <div className="input-group input-group-sm">
          <UseStructureInputAddon value={smarts} selection={selection} setValue={setSmarts} />
          <input
            type="text"
            className="form-control form-control-sm"
            id="smartsStructureInput"
            aria-describedby="smartsStructureInput"
            value={smarts}
            onChange={(e) => setSmarts(e.currentTarget.value)}
          />
        </div>
        <small id="smartsStructureInputHelp" className="form-text text-muted">
          Use a SMILES or SMARTS query
        </small>
      </div>
      <div className="text-end">
        <ButtonWithUpload loading={loading} disabled={!smarts} text="Compute Substructure Matches" />
      </div>
    </FormWrapper>
  );
};
