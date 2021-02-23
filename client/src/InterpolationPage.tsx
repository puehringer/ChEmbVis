import React from "react";
import { interpolateStructures } from "./utils/api";
import { LoadingPage } from "./components/LoadingPage";
import { StructureCardGrid } from "./components/StructureCardGrid";
import { StructureImage } from "./components/StructureImage";
import { ICollection, IInterpolatedParticle } from "./interfaces";
import { HorizontalCollapse } from "./components/HorizontalCollapse";

export function InterpolationPage({
  structures,
  setStructures,
  collection,
  setCollection,
}: {
  structures?: string[];
  setStructures?(structures: string[]): void;
  collection?: ICollection<IInterpolatedParticle>;
  setCollection(collection: ICollection<IInterpolatedParticle>): void;
}) {
  const [optionsCollapsed, setOptionsCollapsed] = React.useState<boolean>(false);
  const [input, setInput] = React.useState<string>(structures?.join("\n") || "");
  const [loading, setLoading] = React.useState<boolean>(false);

  React.useEffect(() => {
    if (structures && structures.length > 0) {
      setLoading(true);
      interpolateStructures(structures, 20).then((data) => {
        setCollection({
          data,
          name: "Interpolated",
        });
        setLoading(false);
      });
    } else {
      setCollection({
        data: [],
        name: "Interpolated",
      });
      setLoading(false);
    }
  }, [structures, setCollection]);

  return (
    <>
      <HorizontalCollapse
        label="Options"
        position="left"
        size="col-md-2"
        collapsed={optionsCollapsed}
        setCollapsed={setOptionsCollapsed}
      >
        <p className="lead">Use the continuous chemical space to interpolate between multiple structures.</p>
        <form
          onSubmit={(e) => {
            e.preventDefault();
            e.stopPropagation();
            setStructures?.(input.split("\n").filter(Boolean));
          }}
        >
          <div className="form-group">
            <label htmlFor="inpolationStructures">Interpolate between:</label>
            <textarea
              className="form-control form-control-sm"
              id="inpolationStructures"
              rows={3}
              onChange={(e) => setInput(e.currentTarget.value)}
              value={input}
            />
          </div>
          <div className="text-right">
            <button className="btn btn-primary" type="submit" disabled={loading}>
              {loading ? (
                <>
                  <span className="spinner-grow spinner-grow-sm" role="status" aria-hidden="true"></span> Loading...
                </>
              ) : (
                <>Go</>
              )}
            </button>
          </div>
          <div className="text-center">
            {structures && structures.length > 0
              ? structures
                  .map<React.ReactNode>((structure, i) => (
                    <StructureImage key={i} structure={structure} width="70px" height="70px" />
                  ))
                  .reduce((prev, curr) => [
                    prev,
                    <i className="fas  fa-fw fa-long-arrow-alt-right ml-2 mr-2"></i>,
                    curr,
                  ])
              : null}
          </div>
        </form>
      </HorizontalCollapse>
      <div
        // className="col-md-10"
        style={{
          display: "flex",
          flex: 1,
          flexDirection: "column",
          overflow: "auto",
          marginLeft: 33,
        }}
      >
        <LoadingPage loading={loading} fallback="Please select structures for interpolation">
          {collection && collection.data.length > 0 ? (
            <StructureCardGrid
              collections={[collection]}
              tableClass="main-ranking"
              structureCardProps={(structure) => ({
                className: structure.scaffold ? "border-primary" : "",
              })}
            />
          ) : null}
        </LoadingPage>
      </div>
    </>
  );
}
