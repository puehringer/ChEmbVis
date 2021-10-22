import * as React from "react";
import { ICollection, IParticle, IParticleSelection } from "../../interfaces";
import { ButtonWithUpload } from "../ButtonWithUpload";
import { getNeighborSamples, getTanimotoSimilarity } from "../../utils/api";
import { FormWrapper } from "./FormWrapper";
import { UseStructureInputAddon } from "../UseStructureInputAddon";
import { Modal } from "react-bootstrap";
import { StructureImage } from "../StructureImage";
import { useNameInput } from "../../utils/hooks";

export const NeighborSamplingForm = ({
  visible,
  setVisible,
  collections,
  setCollections,
  selection,
  loading,
  setLoading,
}: {
  visible: IParticle[] | null;
  setVisible(collection: IParticle[] | null): void;
  collections: ICollection[];
  setCollections(collections: ICollection[]): void;
  selection?: IParticleSelection;
  loading: boolean;
  setLoading(loading: boolean): void;
}) => {
  const [smiles, setSmiles] = React.useState<string>("");
  const [method, setMethod] = React.useState<string>("chembl_pca");
  const [nr, setNr] = React.useState<number>(5);
  const [scale, setScale] = React.useState<number>(1);
  const [name, setName, nameInput] = useNameInput("neighborhoodNameInput", `${method} ${nr}*${scale} Neighborhood of ${smiles}`);

  const computeNeighborSamples = async () => {
    const serverCollection = await getNeighborSamples(smiles, nr, method, scale);
    const similarity = await getTanimotoSimilarity(
      serverCollection.data.map((p) => p.structure),
      smiles,
      "ecfp4"
    );

    const data = serverCollection.data.map((p) => ({
      ...p,
      properties: { ...(p.properties || {}), neighborhoodSimilarity: similarity.tanimoto[p.structure] },
    }));
    setCollections([
      ...collections,
      {
        ...serverCollection,
        name,
        data,
        type: "neighborhoodSampling",
      },
    ]);
    setName("");
    setVisible(data);
  };

  return (
    <>
      <Modal show={Boolean(visible)} onHide={() => setVisible(null)} size="xl">
        <Modal.Header closeButton>
          <Modal.Title>Neighborhood Sampling</Modal.Title>
        </Modal.Header>
        <Modal.Body>
          {visible ? (
            <>
              <div
                style={{
                  display: "grid",
                  gridTemplateColumns: `repeat(${Math.sqrt(visible.length)}, 1fr)`,
                }}
              >
                {visible.map((sample, i) => (
                  <StructureImage
                    key={i}
                    structure={sample.structure}
                    align={visible[Math.floor(visible.length / 2)].structure}
                    title={`Similarity to reference: ${sample.properties?.["neighborhoodSimilarity"]}`}
                    style={{
                      width: "100%",
                      opacity: sample.properties?.["neighborhoodSimilarity"] as number,
                      border: `1px solid ${Math.floor(visible.length / 2) === i ? "gold" : "rgba(0, 0, 0, 0.05)"}`,
                    }}
                  />
                ))}
              </div>
            </>
          ) : null}
        </Modal.Body>
      </Modal>
      <FormWrapper
        title="Neighborhood Sampling"
        loading={loading}
        setLoading={setLoading}
        onSubmit={computeNeighborSamples}
      >
        {nameInput}
        <div className="mb-3">
          <label htmlFor="neighborhoodStructureInput">Reference structure</label>
          <div className="input-group input-group-sm">
            <UseStructureInputAddon value={smiles} selection={selection} setValue={setSmiles} />
            <input
              type="text"
              className="form-control form-control-sm"
              id="neighborhoodStructureInput"
              aria-describedby="neighborhoodStructureInput"
              value={smiles}
              onChange={(e) => setSmiles(e.currentTarget.value)}
            />
          </div>
        </div>
        <div className="mb-3">
          <label htmlFor="neighborhoodSamplesInput" className="form-label">Nr. of Samples: {nr}</label>
          <input
            type="range"
            className="form-range"
            id="neighborhoodSamplesInput"
            value={nr}
            onChange={(e) => setNr(e.currentTarget.valueAsNumber)}
            min={3}
            max={15}
            step={2}
          />
        </div>
        <div className="mb-3">
          <label htmlFor="neighborhoodScaleInput" className="form-label">Scale: {scale}</label>
          <input
            type="range"
            className="form-range"
            id="neighborhoodScaleInput"
            value={scale}
            onChange={(e) => setScale(e.currentTarget.valueAsNumber)}
            min={0.1}
            max={5}
            step={0.1}
          />
        </div>
        <div className="mb-3">
          <label htmlFor="neighborhoodSamplesMethod">Method</label>
          <div className="input-group input-group-sm">
            <select
              className="form-control"
              id="neighborhoodSamplesMethod"
              value={method}
              onChange={(e) => {
                setMethod(e.currentTarget.value);
              }}
            >
              <option value="chembl_pca">PCA of ChEMBL</option>
              <option value="random_orthogonal">Random Orthogonal</option>
            </select>
          </div>
        </div>
        <div className="text-end">
          <ButtonWithUpload loading={loading} disabled={!smiles} text="Compute Neighbors" />
        </div>
      </FormWrapper>
    </>
  );
};
