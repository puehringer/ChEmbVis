import * as React from "react";
import { Modal, Button } from "react-bootstrap";
import { StructureCard } from "./components/StructureCard";
import { StructureCardGrid } from "./components/StructureCardGrid";
import { StructureImage } from "./components/StructureImage";
import { ICollection, IParticle, IParticleSelection } from "./interfaces";
import groupByLodash from "lodash.groupby";

export function getClustersFromParticle(particle: IParticle): string[] {
  return particle?.properties?.clusters?.split(";").filter(Boolean) || [];
}

export const ClusterSidePanel = React.memo(
  ({
    clusters,
    setClusters,
    selected,
    setSelected,
  }: {
    clusters: { [key: string]: IParticle[] };
    setClusters(name: string, clusters: IParticle[] | null): void;
    selected: IParticleSelection;
    setSelected(selected: IParticleSelection): void;
  }) => {
    const selectedParticles = React.useMemo<IParticle[]>(
      () => (selected ? ([] as IParticle[]).concat(...Object.values(selected)) : []),
      [selected]
    );
    const [name, setName] = React.useState<string>("");
    const [wasValidated, setWasValidated] = React.useState<boolean>(false);
    const [showStructures, setShowStructures] = React.useState<
      | {
          type: "selection";
        }
      | {
          type: "cluster";
          cluster: string;
        }
      | null
    >(null);
    const [rankingSelection, setRankingSelection] = React.useState<IParticleSelection>(null);

    const clusterNames = Object.keys(clusters);
    const clusterNameTaken = clusterNames.includes(name);

    const structuresToShow = React.useMemo<ICollection[] | null>(() => {
      if (showStructures) {
        const filteredParticles =
          showStructures.type === "selection"
            ? selectedParticles
            : showStructures.type === "cluster"
            ? clusters[showStructures.cluster]
            : null;
        return filteredParticles
          ? Object.entries(groupByLodash(filteredParticles, "collection")).map(([name, data]) => ({ name, data }))
          : null;
      }
      return null;
    }, [showStructures, clusters, selectedParticles]);

    return (
      <>
        <Modal
          show={Boolean(showStructures && structuresToShow)}
          onHide={() => setShowStructures(null)}
          size="xl"
          dialogClassName="modal-full-width"
        >
          {showStructures && structuresToShow ? (
            <>
              <Modal.Header closeButton>
                <Modal.Title>View {structuresToShow.reduce((a, b) => a + b.data.length, 0)} entries</Modal.Title>
              </Modal.Header>
              <Modal.Body
                style={{
                  display: "flex",
                  flexDirection: "column",
                }}
              >
                <StructureCardGrid collections={structuresToShow} setSelected={(s) => setRankingSelection(s)} />
              </Modal.Body>
              <Modal.Footer>
                <Button
                  variant="primary"
                  onClick={() => {
                    setShowStructures(null);
                    setSelected(rankingSelection);
                  }}
                  disabled={!rankingSelection}
                >
                  Save as new selection
                </Button>
                <Button variant="secondary" onClick={() => setShowStructures(null)}>
                  Close
                </Button>
              </Modal.Footer>
            </>
          ) : null}
        </Modal>
        {selected ? (
          <>
            <p className="lead">Selection</p>
            <p>
              {selectedParticles.length} Selected{" "}
              <button
                type="button"
                className="btn btn-sm btn-light float-right"
                title="Show in table view"
                onClick={() => setShowStructures({ type: "selection" })}
              >
                <i className="fas fa-fw fa-table" />
              </button>
            </p>
          </>
        ) : null}
        {clusterNames.length > 0 ? (
          <>
            <p className="lead">Clusters</p>
            {clusterNames.map((cluster) => (
              <div key={cluster} className="clearfix mb-1">
                {cluster} ({clusters[cluster].length})
                <div className="float-right">
                  <div className="btn-group btn-group-sm ml-2 mr-2" role="group">
                    <button
                      type="button"
                      className="btn btn-light"
                      title="Show cluster in table view"
                      onClick={() => {
                        setShowStructures({ type: "cluster", cluster });
                      }}
                    >
                      <i className="fas fa-fw fa-table" />
                    </button>
                    <button
                      type="button"
                      className="btn btn-danger"
                      title="Delete cluster"
                      onClick={() => {
                        setClusters(cluster, null);
                      }}
                    >
                      <i className="fas fa-fw fa-times" />
                    </button>
                  </div>
                </div>
              </div>
            ))}
          </>
        ) : (
          <p className="lead">No clusters created</p>
        )}
        {selectedParticles.length > 0 ? (
          <>
            <p className="lead">Create new cluster</p>
            {selectedParticles.length === 1 ? (
              <StructureCard structure={selectedParticles[0]} showProperties={false} />
            ) : (
              <StructureImage style={{ maxWidth: "100%" }} structure={selectedParticles.map((s) => s.structure)} />
            )}
            <form
              className={`${wasValidated ? "was-validated" : ""}`}
              onSubmit={(e) => {
                e.preventDefault();
                e.stopPropagation();

                if (clusterNameTaken) {
                  setWasValidated(true);
                  return;
                }
                setClusters(name, selectedParticles);
                setWasValidated(false);
                setName("");
              }}
            >
              <div className="form-group">
                <label htmlFor="clusterNameInput">Cluster name</label>
                <input
                  type="text"
                  className={`form-control form-control-sm ${clusterNameTaken ? "is-invalid" : ""}`}
                  id="clusterNameInput"
                  value={name}
                  required={true}
                  onChange={(e) => setName(e.currentTarget.value)}
                />
                {clusterNameTaken ? (
                  <div className="invalid-feedback">Cluster with this name already exists.</div>
                ) : null}
              </div>

              <div className="text-right">
                <button type="submit" className="btn btn-primary" disabled={clusterNameTaken}>
                  Create cluster
                </button>
              </div>
            </form>
          </>
        ) : (
          <p className="lead">Select to create a new cluster</p>
        )}
      </>
    );
  }
);
