import * as React from "react";
import { ICollection, IParticle, IParticleSelection } from "../../interfaces";
import { ButtonWithUpload } from "../ButtonWithUpload";
import { getNeighborSamples, getTanimotoSimilarity } from "../../utils/api";
import { FormWrapper } from "./FormWrapper";
import { UseStructureInputAddon } from "../UseStructureInputAddon";
import { Modal } from "react-bootstrap";
import { StructureImage } from "../StructureImage";

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
  const [nr, setNr] = React.useState<number>(5);

  const computeNeighborSamples = async () => {
    const particles = await getNeighborSamples(smiles, nr);
    const similarity = await getTanimotoSimilarity(
      particles.map((p) => p.structure),
      smiles,
      "morgan"
    );

    const data = particles.map((p) => ({
      ...p,
      properties: { ...(p.properties || {}), neighborhoodSimilarity: similarity.tanimoto[p.structure] },
    }));
    setCollections([
      ...collections,
      {
        name: `${smiles} Neighborhood`,
        data,
        type: "neighborhoodSampling",
      },
    ]);
    setVisible(data);
  };

  return (
    <>
      <Modal show={visible} onHide={() => setVisible(null)} size="xl">
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
        <div className="form-group">
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
          <div className="form-group">
            <label htmlFor="neighborhoodSamplesInput">Nr. of Samples: {nr}</label>
            <input
              type="range"
              className="form-control-range"
              id="neighborhoodSamplesInput"
              value={nr}
              onChange={(e) => setNr(e.currentTarget.valueAsNumber)}
              min={3}
              max={15}
              step={2}
            />
          </div>
        </div>
        <div className="text-right">
          <ButtonWithUpload loading={loading} disabled={!smiles} text="Compute Similarity" />
        </div>
      </FormWrapper>
    </>
  );
};
