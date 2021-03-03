import React from "react";
import { LoadingPage } from "./components/LoadingPage";
import { StructureCardGrid } from "./components/StructureCardGrid";
import { ICollection, IInterpolatedParticle } from "./interfaces";
import { HorizontalCollapse } from "./components/HorizontalCollapse";
import { InterpolationForm } from "./components/form/InterpolationForm";

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
  const [loading, setLoading] = React.useState<boolean>(false);

  return (
    <>
      <HorizontalCollapse
        label="Options"
        position="left"
        size="col-md-2"
        collapsed={optionsCollapsed}
        setCollapsed={setOptionsCollapsed}
      >
        <InterpolationForm open={true} setCollection={setCollection} loading={loading} setLoading={setLoading} setStructures={setStructures} structures={structures} />
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
