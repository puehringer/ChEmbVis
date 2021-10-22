import * as React from "react";
import { ICollection, IParticle, IParticleSelection } from "../../interfaces";
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
  const [similaryFP, setSimilarityFP] = React.useState<string>("ecfp4");

  const computeTanimotoSimilarity = () => {

    if(similaryFP === "embeddings") {
      setCollections(collections.map((c) => {
        const reference: IParticle = Object.values(selection || {})?.[0]?.[0];
        if(!reference || !reference.embedding || Object.entries(reference.embedding).length === 0) {
          throw Error('This type of similarity requires an active selection (with precomputed embeddings)');
        }
        const availableEmbeddings = Object.keys(reference.embedding);
          
        return {
          ...c,
          data: c.data.map((p) => ({
            ...p,
            properties: {
              ...(p.properties || {}),
              ...(availableEmbeddings.reduce((acc, cur) => ({...acc, [`${cur} Distance ${similarityRef}`]: p.embedding?.[cur]?.reduce((acc, x, i) => {
                const y = reference.embedding![cur][i];
                // TODO: Implement other distance metrics (VAE requires different one...)
                return acc - Math.pow(x - y, 2);
              }, 0)}), {})),
            },
          })),
          selection: null,
        };
      }));
      return;
    }

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
                [`${similaryFP} Distance ${similarityRef}`]: results[i].tanimoto[p.structure],
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
      <div className="mb-3">
        <label htmlFor="similarityRefStructureInput">Reference structure</label>
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
      <div className="mb-3">
        <label htmlFor="similarityFingerprintInput">Fingerprint</label>
        <select
          className="form-control form-control-sm"
          id="similarityFingerprintInput"
          value={similaryFP}
          onChange={(e) => setSimilarityFP(e.currentTarget.value)}
        >
          <option value="ecfp2">ECFP2 (Morgan Fingerprint)</option>
          <option value="ecfp4">ECFP4 (Morgan Fingerprint)</option>
          <option value="cddd">Descriptors from CDDD (ChEMBL only)</option>
          <option value="embeddings">Embeddings</option>
        </select>
      </div>
      <div className="text-end">
        <ButtonWithUpload loading={loading} disabled={!similarityRef} text="Compute Similarity" />
      </div>
    </FormWrapper>
  );
};
