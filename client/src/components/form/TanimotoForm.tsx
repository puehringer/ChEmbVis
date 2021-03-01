import * as React from "react";
import { ICollection, IParticleSelection } from "../../interfaces";
import { ButtonWithUpload } from "../ButtonWithUpload";
import { getTanimotoSimilarity } from "../../utils/api";
import { FormWrapper } from "./FormWrapper";
import { UseStructureInputAddon } from "../UseStructureInputAddon";

export const TanimotoForm = ({
  collections,
  setCollections,
  selection,
  loading,
  setLoading,
}: {
  collections: ICollection[];
  setCollections(collections: ICollection[]): void;
  selection?: IParticleSelection;
  loading: boolean;
  setLoading(loading: boolean): void;
}) => {
  const [similarityRef, setSimilarityRef] = React.useState<string>("");
  const [similaryFP, setSimilarityFP] = React.useState<string>("morgan");

  const computeTanimotoSimilarity = () => {
    return Promise.all(
      collections.map(async ({ data }) =>
        getTanimotoSimilarity(
          data.map((s) => s.structure),
          similarityRef,
          similaryFP
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
                [`Tanimoto ${similaryFP} ${similarityRef}`]: results[i].tanimoto[p.structure],
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
      title="Tanimoto Similarity"
      loading={loading}
      setLoading={setLoading}
      onSubmit={computeTanimotoSimilarity}
    >
      <div className="form-group">
        <label htmlFor="similarityRefStructureInput">Reference structure:</label>
        <div className="input-group input-group-sm">
          <UseStructureInputAddon value={similarityRef} selection={selection} setValue={setSimilarityRef} />
          <input
            type="text"
            className="form-control form-control-sm"
            id="similarityRefStructureInput"
            aria-describedby="similarityRefStructureInput"
            value={similarityRef}
            onChange={(e) => setSimilarityRef(e.currentTarget.value)}
          />
        </div>
      </div>
      <div className="form-group">
        <label htmlFor="similarityFingerprintInput">Fingerprint:</label>
        <select
          className="form-control form-control-sm"
          id="similarityFingerprintInput"
          value={similaryFP}
          onChange={(e) => setSimilarityFP(e.currentTarget.value)}
        >
          <option value="morgan">ECFP4 (Morgan Fingerprint)</option>
          <option value="daylight">Daylight (RDKit Fingerprint)</option>
        </select>
      </div>
      <div className="text-right">
        <ButtonWithUpload loading={loading} disabled={!similarityRef} text="Compute Similarity" />
      </div>
    </FormWrapper>
  );
};
